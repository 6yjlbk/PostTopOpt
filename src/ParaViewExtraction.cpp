#include "ParaViewExtraction.h"
#include <cmath>
#include <stdexcept>
#include <mutex>
#ifdef POSTTOPOPT_WITH_PARAVIEW
#include <vtkAMRDualContour.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCommand.h>
#include <vtkCompositeDataIterator.h>
#include <vtkDataObject.h>
#include <vtkDoubleArray.h>
#include <vtkDummyController.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNew.h>
#include <vtkNonOverlappingAMR.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkTriangleFilter.h>
#include <vtkUniformGrid.h>
#include <vtkVersionMacros.h>
#endif
namespace geometry {
bool paraViewAvailable() noexcept {
#ifdef POSTTOPOPT_WITH_PARAVIEW
    return true;
#else
    return false;
#endif
}
Mesh paraViewDualContour(const CellScalarVolume& volume,const ParaViewContourOptions& options) {
#ifndef POSTTOPOPT_WITH_PARAVIEW
    (void)volume; (void)options;
    throw std::runtime_error("ParaView backend disabled; configure POSTTOPOPT_WITH_PARAVIEW=ON");
#else
    volume.validate();
    for(int n:volume.grid.cells) if(n<3) throw std::invalid_argument("ParaView AMR requires at least 3 cells per axis");
    if(!std::isfinite(options.isovalue)) throw std::invalid_argument("Invalid ParaView isovalue");
    for(double rho:volume.values) if(rho<0 || rho>1) throw std::invalid_argument("ParaView volume fractions must be in [0,1]");
    if(options.isovalue<0 || options.isovalue>1) throw std::invalid_argument("Volume fraction threshold must be in [0,1]");
    // ParaView's AMR helper dereferences the global controller even in serial runs.
    // Restore the host controller on all exits; serialize our adapter calls.
    static std::mutex controllerMutex;
    const std::lock_guard<std::mutex> lock(controllerMutex);
    struct ControllerScope {
        vtkSmartPointer<vtkMultiProcessController> previous;
        vtkNew<vtkDummyController> serial;
        ControllerScope():previous(vtkMultiProcessController::GetGlobalController()) {
            vtkMultiProcessController::SetGlobalController(serial);
        }
        ~ControllerScope(){vtkMultiProcessController::SetGlobalController(previous);}
    } controllerScope;
    const auto& g=volume.grid;
    vtkNew<vtkUniformGrid> block;
    block->SetOrigin(g.origin.data()); block->SetSpacing(g.spacing.data());
    block->SetDimensions(g.cells[0]+1,g.cells[1]+1,g.cells[2]+1);
    vtkNew<vtkDoubleArray> values;
    values->SetName("density"); values->SetNumberOfComponents(1);
    values->SetNumberOfTuples(block->GetNumberOfCells());
    // VTK has x-fastest indexing, PostTopOpt has z-fastest indexing.
    for(int z=0;z<g.cells[2];++z) for(int y=0;y<g.cells[1];++y) for(int x=0;x<g.cells[0];++x) {
        vtkIdType id=(static_cast<vtkIdType>(z)*g.cells[1]+y)*g.cells[0]+x;
        values->SetTuple1(id,volume.values[cellIndex(g,x,y,z)]);
    }
    block->GetCellData()->SetScalars(values);
    vtkNew<vtkNonOverlappingAMR> amr;
#if VTK_VERSION_NUMBER >= VTK_VERSION_CHECK(9,6,0)
    amr->Initialize(std::vector<unsigned int>{1});
#else
    int blocks[1]={1}; amr->Initialize(1,blocks);
#endif
    amr->SetDataSet(0,0,block);
    class Errors final:public vtkCommand {
    public:
        static Errors* New(){return new Errors;}
        bool failed=false;
        void Execute(vtkObject*,unsigned long,void*) override {failed=true;}
    };
    vtkNew<Errors> errors;
    vtkNew<vtkAMRDualContour> contour;
    contour->AddObserver(vtkCommand::ErrorEvent,errors);
    contour->SetInputDataObject(amr);
    contour->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"density");
    contour->SetIsoValue(options.isovalue);
    contour->SetEnableCapping(options.enableCapping);
    contour->SetEnableMergePoints(options.mergePoints);
    vtkNew<vtkDummyController> controller;
    contour->SetController(controller);
    contour->SetEnableMultiProcessCommunication(1); // Serial controller; helper requires a non-null communicator.
    contour->Update();
    if(errors->failed) throw std::runtime_error("vtkAMRDualContour reported an error");
    auto* output=vtkMultiBlockDataSet::SafeDownCast(contour->GetOutputDataObject(0));
    if(!output) throw std::runtime_error("Unexpected ParaView output type");
    Mesh mesh;
    auto iterator=vtkSmartPointer<vtkCompositeDataIterator>::Take(output->NewIterator());
    iterator->SkipEmptyNodesOn();
    for(iterator->InitTraversal();!iterator->IsDoneWithTraversal();iterator->GoToNextItem()) {
        auto* poly=vtkPolyData::SafeDownCast(iterator->GetCurrentDataObject());
        if(!poly) continue;
        vtkNew<vtkTriangleFilter> triangles;
        triangles->AddObserver(vtkCommand::ErrorEvent,errors);
        triangles->SetInputData(poly); triangles->PassLinesOff(); triangles->PassVertsOff(); triangles->Update();
        if(errors->failed) throw std::runtime_error("ParaView output triangulation failed");
        auto* data=triangles->GetOutput();
        auto offset=mesh.vertices.size();
        for(vtkIdType i=0;i<data->GetNumberOfPoints();++i) {
            double p[3]; data->GetPoint(i,p); Vec3 v(p[0],p[1],p[2]);
            if(!v.allFinite()) throw std::runtime_error("Non-finite ParaView vertex");
            mesh.vertices.push_back(v);
        }
        vtkIdType count; const vtkIdType* ids=nullptr;
        auto* cells=data->GetPolys(); cells->InitTraversal();
        while(cells->GetNextCell(count,ids)) {
            if(count!=3) throw std::runtime_error("Unexpected polygon after triangulation");
            for(int j=0;j<3;++j) if(ids[j]<0 || ids[j]>=data->GetNumberOfPoints()) throw std::runtime_error("Invalid ParaView index");
            mesh.triangles.push_back({offset+static_cast<std::size_t>(ids[0]),offset+static_cast<std::size_t>(ids[1]),offset+static_cast<std::size_t>(ids[2])});
        }
    }
    return mesh;
#endif
}
}




