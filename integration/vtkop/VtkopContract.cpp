#include "VtkopContract.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <vector>
#include <vtkAMRDualContour.h>
#include <vtkBoxClipDataSet.h>
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellDataToPointData.h>
#include <vtkCellType.h>
#include <vtkCommand.h>
#include <vtkCompositeDataIterator.h>
#include <vtkDataArray.h>
#include <vtkDataObject.h>
#include <vtkDataSetAttributes.h>
#include <vtkDoubleArray.h>
#include <vtkDummyController.h>
#include <vtkExtractCells.h>
#include <vtkIdList.h>
#include <vtkImageData.h>
#include <vtkMergePoints.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkMultiProcessController.h>
#include <vtkNew.h>
#include <vtkNonOverlappingAMR.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkProbeFilter.h>
#include <vtkSmartPointer.h>
#include <vtkStaticCellLocator.h>
#include <vtkTetra.h>
#include <vtkTriangleFilter.h>
#include <vtkUniformGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVersionMacros.h>
// Reference implementation for the proposed vtkop contract.
// Requires a compatible ParaView SDK, not standalone VTK.
// Integrate these wrappers into the real vtkop;
// no application pipeline here.
struct AmrObject {
    vtkSmartPointer<vtkNonOverlappingAMR> data =
    vtkSmartPointer<vtkNonOverlappingAMR>::New();
};

struct ContourObject {
    vtkSmartPointer<vtkAMRDualContour> filter = vtkSmartPointer<vtkAMRDualContour>::New();
    // РЎРѕР·РґР°РЅРёРµ РѕР±СЉРµРєС‚РІ VTK-С„РёР»СЊС‚СЂР° РґР»СЏ DC
    Amr input;
    bool updated = false;
    // РЎР»СѓР¶РµР±РЅС‹Р№ С„Р»Р°Рі РІС‹РїРѕР»РЅРµРЅРёСЏ СЂР°Р±РѕС‚С‹ Р°Р»РіРѕСЂРёС‚РјР°
};

struct PolyDataObject {
    vtkSmartPointer<vtkPolyData> data;
    // РЎРѕР·РґР°РЅРёРµ РѕР±СЉРµРєС‚РІ РґР»СЏ С…СЂР°РЅРµРЅРёСЏ РґР°РЅРЅС‹С… VTK
    bool triangulated = false;
    // РЎР»СѓР¶РµР±РЅС‹Р№ С„Р»Р°Рі РІС‹РїРѕР»РЅРµРЅРёСЏ СЂР°Р±РѕС‚С‹ Р°Р»РіРѕСЂРёС‚РјР°
};

// Р‘РµР·С‹РјСЏРЅРЅРѕРµ РїСЂРѕСЃС‚СЂР°РЅСЃС‚РІРѕ РёРјС‘РЅ, РІ РєРѕС‚РѕСЂРѕРј РЅР°С…РѕРґРёС‚СЃСЏ РІСЃРїРѕРјРѕРіР°С‚РµР»СЊРЅС‹Р№ РєРѕРґ (РїСЂРѕРІРµСЂРєР°, СЂР°Р±РѕС‚Р° СЃ РѕС€РёР±РєР°РјРё, СЂР°Р±РѕС‚Р° СЃ С„РѕСЂРјР°С‚Р°РјРё РґР°РЅРЅС‹С…)
namespace {
    template<class T> T& checked(const std::shared_ptr<T>& object) {
        // РџСЂРѕРІРµСЂРєР°, С‡С‚Рѕ СѓРєР°Р·Р°С‚РµР»СЊ РЅР° РѕР±СЉРµРєС‚ РЅРµ РїСѓСЃС‚
        if (!object) throw std::invalid_argument("Null vtkop handle");
        return *object;
    }
    // РџРµСЂРµРІРѕРґ СЂР°Р·РјРµСЂРѕРІ РґР»СЏ VTK
    vtkIdType vtkSize(std::size_t value) {
        if (static_cast<std::uintmax_t>(value) >
        static_cast<std::uintmax_t>((std::numeric_limits<vtkIdType>::max)()))
        throw std::runtime_error("Size exceeds vtkIdType range");
        return static_cast<vtkIdType>(value);
    }
    // РћР±СЂР°С‚РЅРѕРµ РїСЂРµРѕР±СЂР°Р·РѕРІР°РЅРёРµ СЂР°Р·РјРµСЂРѕРІ
    std::size_t nativeSize(vtkIdType value) {
        if (value < 0 || static_cast<std::uintmax_t>(value) >
        static_cast<std::uintmax_t>((std::numeric_limits<std::size_t>::max)()))
        throw std::runtime_error("Invalid VTK size or index");
        return static_cast<std::size_t>(value);
    }
    class Errors final : public vtkCommand { // Р—Р°РїРѕРјРёРЅР°РµС‚ РѕС€РёР±РєСѓ, С‡С‚РѕР±С‹ РїРѕСЃР»Рµ РІС‹Р·РѕРІР° С„РёР»СЊС‚СЂР° СЃРѕРѕР±С‰РёС‚СЊ Рѕ РЅРµР№
        public:
        static Errors* New() {
            return new Errors;
        }
        bool failed = false;
        void Execute(vtkObject*, unsigned long, void*) override { failed = true;
        }
    };
    struct ObserveErrors { // РђРІС‚РѕРјР°С‚РёС‡РµСЃРєРё РїРѕРґРєР»СЋС‡Р°РµС‚ Рё РѕС‚РєР»СЋС‡Р°РµС‚ РѕР±СЂР°Р±РѕС‚С‡РёРє РѕС€РёР±РѕРє
        vtkObject* object;
        unsigned long tag;
        ObserveErrors(vtkObject* o, vtkCommand* errors)
        : object(o), tag(o->AddObserver(vtkCommand::ErrorEvent, errors)) {
        }
        ~ObserveErrors() {
            object->RemoveObserver(tag);
        }
    };
    struct ControllerScope { // Р’СЂРµРјРµРЅРЅРѕ РјРµРЅСЏРµС‚ РіР»РѕР±Р°Р»СЊРЅС‹Р№ РєРѕРЅС‚СЂРѕР»Р»РµСЂ Рё РіР°СЂР°РЅС‚РёСЂСѓРµС‚ РІРѕСЃСЃС‚Р°РЅРѕРІР»РµРЅРёРµ РїСЂРµРґС‹РґСѓС‰РµРіРѕ
        vtkSmartPointer<vtkMultiProcessController> previous;
        vtkNew<vtkDummyController> serial;
        ControllerScope() : previous(vtkMultiProcessController::GetGlobalController()) {
            vtkMultiProcessController::SetGlobalController(serial);
        }
        ~ControllerScope() {
            vtkMultiProcessController::SetGlobalController(previous);
        }
    };
    constexpr const char* Density="density";
    // vtkBoxClipDataSet creates its own vtkPoints as float. The insertion
    // locator is initialized before any output points are inserted: select double
    // there, retaining vtkMergePoints connectivity and duplicate-point handling.
    class DoublePrecisionMergePoints final : public vtkMergePoints {
        public:
        static DoublePrecisionMergePoints* New() {
            return new DoublePrecisionMergePoints;
        }
        int InitPointInsertion(vtkPoints* points,const double bounds[6]) override {
            points->SetDataTypeToDouble();
            return this->vtkMergePoints::InitPointInsertion(points,bounds);
        }
        int InitPointInsertion(vtkPoints* points,const double bounds[6],vtkIdType estimate) override {
            points->SetDataTypeToDouble();
            return this->vtkMergePoints::InitPointInsertion(points,bounds,estimate);
        }
    };
    void require(bool ok,const char* message) {
        if (!ok)throw std::runtime_error(message);
    }
    double tetraVolume(vtkCell* cell) {
        require(cell && cell->GetCellType()==VTK_TETRA,"Expected a linear tetrahedron");
        double p[4][3];
        for (int k=0;k<4;++k)cell->GetPoints()->GetPoint(k,p[k]);
        return std::abs(vtkTetra::ComputeVolume(p[0],p[1],p[2],p[3]));
    }
}

TransferGrid makeTransferGrid(vtkUnstructuredGrid* source,
    const std::array<int,3>& cells) {
    require(source && source->GetNumberOfPoints()>0,"Source has no points");
    double bounds[6];
    source->GetBounds(bounds);
    TransferGrid grid;
    for (int axis = 0; axis < 3; ++axis) {
        const double extent = bounds[2 * axis + 1] - bounds[2 * axis];
        require(std::isfinite(extent) && extent > 0 && cells[axis] > 0 &&
        cells[axis] < (std::numeric_limits<int>::max)(), "Invalid bounds or resolution");
        grid.origin[axis] = bounds[2 * axis];
        grid.spacing[axis] = extent / cells[axis];
        grid.cells[axis] = cells[axis];
    }
    return grid;
}

vtkSmartPointer<vtkNonOverlappingAMR> transferVF(vtkUnstructuredGrid* source,
    const std::string& name,const TransferGrid& grid,TransferMethod method) {
    require(source && source->GetNumberOfCells()>0,"Source must contain volume cells");
    require(method==TransferMethod::CellCenterSampling || method==TransferMethod::NodalInterpolation || method==TransferMethod::VolumeWeighted,"Unknown transfer method");
    vtkIdType count=1;
    for (int d=0;d<3;++d) {
        require(grid.cells[d]>0 && grid.cells[d]<std::numeric_limits<int>::max(),"Invalid target dimensions");
        require(std::isfinite(grid.origin[d]) && std::isfinite(grid.spacing[d]) && grid.spacing[d]>0,"Invalid target geometry");
        require(std::isfinite(grid.origin[d]+grid.cells[d]*grid.spacing[d]),"Target bounds overflow");
        require(count<=std::numeric_limits<vtkIdType>::max()/grid.cells[d],"Target cell count overflow");
        count*=grid.cells[d];
    }
    auto vf=source->GetCellData()->GetArray(name.c_str());
    require(vf && vf->GetNumberOfComponents()==1 && vf->GetNumberOfTuples()==source->GetNumberOfCells(),"Expected scalar cell VF array");
    for (vtkIdType i=0;i<source->GetNumberOfPoints();++i) {
        double p[3];
        source->GetPoint(i,p);
        for (double v:p)require(std::isfinite(v),"Nonfinite source point");
    }
    // Isolate the selected array: source point data must not shadow the cell VF.
    vtkNew<vtkUnstructuredGrid> clean;
    clean->CopyStructure(source);
    for (vtkIdType i=0;i<vf->GetNumberOfTuples();++i) {
        double v=vf->GetComponent(i,0);
        require(std::isfinite(v)&&v>=0&&v<=1,"VF must be finite and in [0,1]");
        auto cell=source->GetCell(i);
        require(cell && cell->GetCellDimension()==3 && cell->IsLinear(),"Only linear 3D source cells supported");
        if (method==TransferMethod::VolumeWeighted) {
            require(cell->GetCellType()==VTK_TETRA || cell->GetCellType()==VTK_VOXEL,"VolumeWeighted currently supports TETRA and VOXEL only");
            if (cell->GetCellType()==VTK_TETRA)require(tetraVolume(cell)>0,"Degenerate tetrahedron");
            else {
                double b[6];
                cell->GetBounds(b);
                require(b[1]>b[0]&&b[3]>b[2]&&b[5]>b[4],"Degenerate voxel");
            }
        }
    }
    clean->GetCellData()->SetScalars(vf);
    vtkNew<vtkUniformGrid> block;
    block->SetOrigin(grid.origin.data());
    block->SetSpacing(grid.spacing.data());
    block->SetDimensions(grid.cells[0]+1,grid.cells[1]+1,grid.cells[2]+1);
    vtkNew<vtkDoubleArray> result;
    result->SetName(Density);
    result->SetNumberOfTuples(count);
    result->FillValue(0);
    vtkNew<Errors> errors;
    if (method!=TransferMethod::VolumeWeighted) {
        vtkNew<vtkImageData> centers;
        centers->SetOrigin(grid.origin[0]+.5*grid.spacing[0],grid.origin[1]+.5*grid.spacing[1],grid.origin[2]+.5*grid.spacing[2]);
        centers->SetSpacing(grid.spacing.data());
        centers->SetDimensions(grid.cells.data());
        vtkNew<vtkCellDataToPointData> nodes;
        vtkNew<vtkProbeFilter> probe;
        probe->AddObserver(vtkCommand::ErrorEvent,errors);
        probe->SetInputData(centers);
        if (method==TransferMethod::NodalInterpolation) {
            nodes->AddObserver(vtkCommand::ErrorEvent,errors);
            nodes->SetInputData(clean);
            nodes->Update();
            require(!errors->failed,"CellDataToPointData failed");
            probe->SetSourceData(nodes->GetOutput());
        }
        else probe->SetSourceData(clean);
        vtkNew<vtkStaticCellLocator> locator;
        probe->SetCellLocatorPrototype(locator);
        probe->Update();
        require(!errors->failed,"ProbeFilter failed");
        auto pd=probe->GetOutput()->GetPointData();
        auto sampled=pd->GetArray(name.c_str());
        auto mask=pd->GetArray("vtkValidPointMask");
        require(sampled&&mask&&sampled->GetNumberOfTuples()==count&&mask->GetNumberOfTuples()==count,"Missing probe output or validity mask");
        for (vtkIdType i = 0; i < count; ++i) {
            if (mask->GetComponent(i, 0) == 0) continue;
            double v=sampled->GetComponent(i,0);
            require(std::isfinite(v), "Nonfinite sampled VF");
            result->SetValue(i,std::clamp(v,0.,1.));
        }
    }
    else {
        const double cellVolume=grid.spacing[0]*grid.spacing[1]*grid.spacing[2];
        require(std::isfinite(cellVolume)&&cellVolume>0,"Target volume overflow/underflow");
        vtkNew<vtkStaticCellLocator> locator;
        locator->SetDataSet(clean);
        locator->BuildLocator();
        vtkNew<vtkIdList> ids;
        vtkNew<vtkExtractCells> extract;
        extract->SetInputData(clean);
        extract->AddObserver(vtkCommand::ErrorEvent,errors);
        vtkNew<DoublePrecisionMergePoints> clipLocator;
        vtkNew<vtkBoxClipDataSet> clip;
        clip->SetLocator(clipLocator);
        clip->SetOrientation(0);
        clip->GenerateClippedOutputOff();
        clip->AddObserver(vtkCommand::ErrorEvent,errors);
        for (int z = 0; z < grid.cells[2]; ++z)
        for (int y = 0; y < grid.cells[1]; ++y)
            for (int x = 0; x < grid.cells[0]; ++x) {
            const int index[3]={x,y,z};
            double b[6];
            for (int d=0;d<3;++d) {
                b[2*d]=grid.origin[d]+index[d]*grid.spacing[d];
                b[2*d+1]=grid.origin[d]+(index[d]+1)*grid.spacing[d];
            }
            ids->Reset();
            locator->FindCellsWithinBounds(b,ids);
            if (ids->GetNumberOfIds()==0)continue;
            extract->SetCellList(ids);
            extract->Modified();
            extract->Update();
            require(!errors->failed,"ExtractCells failed");
            clip->SetInputData(extract->GetOutput());
            clip->SetBoxClip(b[0],b[1],b[2],b[3],b[4],b[5]);
            clip->Modified();
            clip->Update();
            require(!errors->failed,"BoxClipDataSet failed");
            auto pieces=clip->GetOutput();
            require(pieces->GetNumberOfPoints()==0 || pieces->GetPoints()->GetDataType()==VTK_DOUBLE,"Clipping must preserve double precision");
            auto values=pieces->GetCellData()->GetArray(name.c_str());
            require(pieces->GetNumberOfCells()==0||values,"Clipping lost VF");
            double volume=0,mass=0;
            for (vtkIdType k=0;k<pieces->GetNumberOfCells();++k) {
                double v=tetraVolume(pieces->GetCell(k));
                volume+=v;
                mass+=v*values->GetComponent(k,0);
            }
            if (!(std::isfinite(volume)&&std::isfinite(mass)&&volume<=cellVolume*(1+1e-6))) {
                throw std::runtime_error("Invalid intersection volume: overlapping cells or clipping error");
            }
            vtkIdType id=(vtkIdType(z)*grid.cells[1]+y)*grid.cells[0]+x;
            result->SetValue(id,std::clamp(mass/cellVolume,0.,1.));
        }
    }
    block->GetCellData()->SetScalars(result);
    auto amr=vtkSmartPointer<vtkNonOverlappingAMR>::New();
#if VTK_VERSION_NUMBER >= VTK_VERSION_CHECK(9,6,0)
    amr->Initialize(std::vector<unsigned int>{1});
#else
    int blocks[1]={1};
    amr->Initialize(1,blocks);
#endif
    amr->SetDataSet(0,0,block);
    return amr;
}

Amr transferCellDensities(const std::vector<std::array<double, 3>>& points,
const std::vector<std::size_t>& offsets,
const std::vector<std::size_t>& nodes,
const std::vector<double>& values,
const std::array<int, 3>& searchCells,
TransferMethod method) {
    if (points.empty() || values.empty() || offsets.empty() ||
    offsets.size() - 1 != values.size() || offsets.front() != 0 || offsets.back() != nodes.size())
    throw std::invalid_argument("Invalid FE connectivity or VF count");
    vtkNew<vtkPoints> coordinates;
    coordinates->SetDataTypeToDouble();
    coordinates->SetNumberOfPoints(vtkSize(points.size()));
    for (std::size_t index = 0; index < points.size(); ++index)
    coordinates->SetPoint(vtkSize(index), points[index].data());
    vtkNew<vtkUnstructuredGrid> source;
    source->SetPoints(coordinates);
    source->Allocate(vtkSize(values.size()));
    vtkNew<vtkDoubleArray> densities;
    densities->SetName("density");
    densities->SetNumberOfTuples(vtkSize(values.size()));
    for (std::size_t cell = 0; cell < values.size(); ++cell) {
        if (offsets[cell] > offsets[cell + 1] || offsets[cell + 1] > nodes.size())
        throw std::invalid_argument("Invalid FE offsets");
        const auto count = offsets[cell + 1] - offsets[cell];
        int type;
        switch (count) {
            case 4: type = VTK_TETRA;
            break;
            case 6: type = VTK_WEDGE;
            break;
            case 8: type = VTK_HEXAHEDRON;
            break;
            default: throw std::invalid_argument("Expected linear solid tetra, wedge or hexahedron");
        }
        vtkIdType ids[8];
        for (std::size_t local = 0; local < count; ++local) {
            const auto node = nodes[offsets[cell] + local];
            if (node >= points.size()) throw std::invalid_argument("FE node index out of range");
            ids[local] = vtkSize(node);
        }
        source->InsertNextCell(type, static_cast<vtkIdType>(count), ids);
        densities->SetValue(vtkSize(cell), values[cell]);
    }
    source->GetCellData()->SetScalars(densities);
    const auto grid = makeTransferGrid(source, searchCells);
    auto result = std::make_shared<AmrObject>();
    result->data = transferVF(source, "density", grid, method);
    return result;
}

// РЎРѕР·РґР°РЅРёРµ РѕР±СЉРµРєС‚Р° С„РёР»СЊС‚СЂР° DC
Contour createAmrDualContour() {
    return std::make_shared<ContourObject>();
}

// РџРµСЂРµРґР°С‡Р° СЏС‡РµРµРє СЃ РІРёСЂС‚СѓР°Р»РЅРѕР№ РїР»РѕС‚РЅРѕСЃС‚СЊСЋ
void setDensityInput(const Contour& contour, const Amr& input) {
    auto& contourObject = checked(contour);
    auto& amr = checked(input);
    auto* block = amr.data->GetDataSet(0, 0);
    if (!block || !block->GetCellData()->GetScalars())
    throw std::invalid_argument("Input cell densities are missing");
    contourObject.updated = false;
    contourObject.input = input;
    contourObject.filter->SetInputDataObject(amr.data);
    contourObject.filter->SetInputArrayToProcess(0, 0, 0,
    vtkDataObject::FIELD_ASSOCIATION_CELLS, Density);
}

// Р—Р°РґР°РЅРёРµ СѓСЂРѕРІРЅСЏ Iso-РїРѕРІРµСЂС…РЅРѕСЃС‚Рё
void setIsovalue(const Contour& contour, double value) {
    auto& contourObject = checked(contour);
    if (!std::isfinite(value) || value < 0 || value > 1)
    throw std::invalid_argument("Invalid VF iso");
    contourObject.updated = false;
    contourObject.filter->SetIsoValue(value);
}

// Р’РєР»/РІС‹РєР» Р·Р°РєСЂС‹С‚РёРµ РїРѕРІРµСЂС…РЅРѕСЃС‚Рё РЅР° РіСЂР°РЅРёС†Рµ СЂР°СЃС‡РµС‚РЅРѕР№ РѕР±Р»Р°СЃС‚Рё
void setCapping(const Contour& contour, bool enabled) {
    auto& contourObject = checked(contour);
    contourObject.updated = false;
    contourObject.filter->SetEnableCapping(enabled);
}

// Р’РєР»/РІС‹РєР» РѕР±СЉРµРґРёРЅРµРЅРёРµ РѕР±С‰РёС… С‚РѕС‡РµРє РІС‹С…РѕРґРЅРѕР№ СЃРµС‚РєРё
void setMergePoints(const Contour& contour, bool enabled) {
    auto& contourObject = checked(contour);
    contourObject.updated = false;
    contourObject.filter->SetEnableMergePoints(enabled);
}

// Р—Р°РїСѓСЃРє Р°Р»РіРѕСЂРёС‚РјР°
void update(const Contour& contour, double* extractionSeconds) {
    if (extractionSeconds) *extractionSeconds = 0;
    auto& contourObject = checked(contour);
    contourObject.updated = false;
    if (!contourObject.input)
    throw std::runtime_error("DC input or cell densities are missing");
    // Serializes THIS implementation only. The real vtkop must coordinate
    // this scope with other host users of the global VTK controller.
    static std::mutex controllerMutex;
    const std::lock_guard<std::mutex> lock(controllerMutex);
    ControllerScope controller;
    contourObject.filter->SetController(controller.serial);
    contourObject.filter->SetEnableMultiProcessCommunication(1);
    vtkNew<Errors> errors;
    ObserveErrors observer(contourObject.filter, errors);
    const auto start = std::chrono::steady_clock::now();
    contourObject.filter->Update();
    if (extractionSeconds)
        *extractionSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    if (errors->failed) throw std::runtime_error("vtkAMRDualContour failed");
    if (!vtkMultiBlockDataSet::SafeDownCast(contourObject.filter->GetOutputDataObject(0)))
    throw std::runtime_error("Unexpected DC output type");
    contourObject.updated = true;
}

// РљРѕРЅС‚РµР№РЅРµСЂ СЃ СЂРµР·СѓР»СЊС‚Р°С‚Р°РјРё DC (СЃРїРёСЃРѕРє РїРѕР»РёРіРѕРЅР°Р»СЊРЅС‹С… Р±Р»РѕРєРѕРІ СЂРµР·СѓР»СЊС‚Р°С‚Р° СЂР°Р±РѕС‚С‹ Р°Р»РіРѕСЂРёС‚РјР° DC)
std::vector<PolyData> polyDataBlocks(const Contour& contour) {
    auto& contourObject = checked(contour);
    if (!contourObject.updated) throw std::runtime_error("DC must be updated before reading output");
    auto* output = vtkMultiBlockDataSet::SafeDownCast(contourObject.filter->GetOutputDataObject(0));
    if (!output) throw std::runtime_error("Missing DC output");
    std::vector<PolyData> result;
    auto iterator = vtkSmartPointer<vtkCompositeDataIterator>::Take(output->NewIterator());
    iterator->SkipEmptyNodesOn();
    for (iterator->InitTraversal(); !iterator->IsDoneWithTraversal(); iterator->GoToNextItem()) {
        auto* poly = vtkPolyData::SafeDownCast(iterator->GetCurrentDataObject());
        if (!poly) continue;
        auto item = std::make_shared<PolyDataObject>();
        item->data = poly;
        result.push_back(item);
    }
    return result;
}

// РўСЂРёР°РЅРіСѓР»РёСЂРѕРІР°РЅРёРµ СЂРµР·СѓР»СЊС‚Р°С‚РѕРІ СЂР°Р±РѕС‚С‹ Р°Р»РіРѕСЂРёС‚РјР° (РїРѕР»РёРіРѕРЅР°Р»СЊРЅР°СЏ РјРѕРґРµР»СЊ -> С‚СЂРёР°РЅРіСѓР»РёСЂРѕРІР°РЅРёРµ РїСЂРµРґСЃС‚Р°РІР»РµРЅРёРµ)
PolyData triangulate(const PolyData& input) {
    auto& source = checked(input);
    vtkNew<vtkTriangleFilter> filter;
    vtkNew<Errors> errors;
    ObserveErrors observer(filter, errors);
    filter->SetInputData(source.data);
    filter->PassLinesOff();
    filter->PassVertsOff();
    filter->Update();
    if (errors->failed || !filter->GetOutput())
    throw std::runtime_error("VTK triangulation failed");
    auto result = std::make_shared<PolyDataObject>();
    result->data = filter->GetOutput();
    result->triangulated = true;
    return result;
}

std::size_t pointCount(const PolyData& data) {
    return nativeSize(checked(data).data->GetNumberOfPoints());
}

std::size_t triangleCount(const PolyData& data) {
    auto& p = checked(data);
    if (!p.triangulated) throw std::runtime_error("Triangulate before reading faces");
    return nativeSize(p.data->GetNumberOfPolys());
}

std::array<double, 3> point(const PolyData& data, std::size_t index) {
    auto& p = checked(data);
    if (index >= pointCount(data)) throw std::out_of_range("VTK point index");
    std::array<double, 3> result{};
    p.data->GetPoint(vtkSize(index), result.data());
    for (double coordinate : result)
    if (!std::isfinite(coordinate))
    throw std::runtime_error("Non-finite surface coordinate");
    return result;
}

std::array<std::size_t, 3> triangle(const PolyData& data, std::size_t index) {
    auto& p = checked(data);
    if (!p.triangulated) throw std::runtime_error("Triangulate before reading faces");
    if (index >= nativeSize(p.data->GetNumberOfPolys()))
    throw std::out_of_range("VTK triangle index");
    vtkIdType count;
    const vtkIdType* ids;
    p.data->GetPolys()->GetCellAtId(vtkSize(index), count, ids);
    if (count != 3) throw std::runtime_error("Expected a triangle");
    return {nativeSize(ids[0]), nativeSize(ids[1]), nativeSize(ids[2])};
}
