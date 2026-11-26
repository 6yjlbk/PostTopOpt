#include "ParaViewExtraction.h"
#include <cmath>
#include <chrono>
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
#include <vtkQuadricDecimation.h>
#include <vtkSmartPointer.h>
#include <vtkTriangleFilter.h>
#include <vtkUniformGrid.h>
#include <vtkVersionMacros.h>
#include <vtkWindowedSincPolyDataFilter.h>
#endif
namespace geometry {
#ifdef POSTTOPOPT_WITH_PARAVIEW
namespace {
// Manual, compile-time postprocessing settings. Rebuild after changing them.
// Both false retain the original triangulation-only path. ParaView backend only.
constexpr bool EnableSmoothing = false; //Сглаживаение вкл/выкл
constexpr bool EnableDecimation = false; //Облегчение модели за счет уменьшения числа граней вкл/выкл
constexpr int SmoothingIterations = 20;
constexpr double SmoothingPassBand = 0.1; // Smaller means stronger smoothing (0,2].
constexpr bool SmoothingBoundary = false; // Keep open-boundary vertices fixed.
constexpr bool SmoothingFeatureEdges = true; // Classify sharp edges and constrain motion.
constexpr double SmoothingFeatureAngle = 45.0; // Degrees.
constexpr double SmoothingEdgeAngle = 15.0; // Degrees.
constexpr bool SmoothingNormalizeCoordinates = true;
constexpr bool SmoothingNonManifold = false;
constexpr double DecimationTargetReduction = 0.5; // Fraction to REMOVE, not retain.
constexpr bool DecimationVolumePreservation = true; // Not an exact volume guarantee.
constexpr bool DecimationAttributeErrorMetric = false;

static_assert(SmoothingIterations >= 2, "Use at least two smoothing iterations");
static_assert(SmoothingPassBand > 0.0 && SmoothingPassBand <= 2.0, "Invalid passband");
static_assert(SmoothingFeatureAngle >= 0.0 && SmoothingFeatureAngle <= 180.0,
              "Invalid feature angle");
static_assert(SmoothingEdgeAngle >= 0.0 && SmoothingEdgeAngle <= 180.0,
              "Invalid edge angle");
static_assert(DecimationTargetReduction >= 0.0 && DecimationTargetReduction < 1.0,
              "Reduction must be in [0,1)");

class Errors final:public vtkCommand {
public:
    static Errors* New(){return new Errors;}
    bool failed=false;
    void Execute(vtkObject*,unsigned long,void*) override {failed=true;}
};

// Retain each output with a smart pointer after its producing filter is destroyed.
// SetInputData does not modify the original contour output in place.
vtkSmartPointer<vtkPolyData> postprocessPolyData(vtkPolyData* input, Errors* errors) {
    vtkSmartPointer<vtkPolyData> data = input;
    if constexpr(EnableSmoothing) {
        if(data->GetNumberOfPoints() > 0 && data->GetNumberOfCells() > 0) {
            vtkNew<vtkWindowedSincPolyDataFilter> smoothing;
            smoothing->AddObserver(vtkCommand::ErrorEvent,errors);
            smoothing->SetInputData(data);
            smoothing->SetNumberOfIterations(SmoothingIterations);
            smoothing->SetPassBand(SmoothingPassBand);
            smoothing->SetBoundarySmoothing(SmoothingBoundary);
            smoothing->SetFeatureEdgeSmoothing(SmoothingFeatureEdges);
            smoothing->SetFeatureAngle(SmoothingFeatureAngle);
            smoothing->SetEdgeAngle(SmoothingEdgeAngle);
            smoothing->SetNormalizeCoordinates(SmoothingNormalizeCoordinates);
            smoothing->SetNonManifoldSmoothing(SmoothingNonManifold);
            smoothing->GenerateErrorScalarsOff();
            smoothing->GenerateErrorVectorsOff();
            smoothing->Update();
            if(errors->failed) throw std::runtime_error("ParaView output smoothing failed");
            data = smoothing->GetOutput();
        }
    }
    vtkNew<vtkTriangleFilter> triangles;
    triangles->AddObserver(vtkCommand::ErrorEvent,errors);
    triangles->SetInputData(data);
    triangles->PassLinesOff(); triangles->PassVertsOff(); triangles->Update();
    if(errors->failed) throw std::runtime_error("ParaView output triangulation failed");
    data = triangles->GetOutput();
    if constexpr(EnableDecimation) {
        if(data->GetNumberOfPolys() > 0 && DecimationTargetReduction > 0.0) {
            vtkNew<vtkQuadricDecimation> decimation;
            decimation->AddObserver(vtkCommand::ErrorEvent,errors);
            decimation->SetInputData(data);
            decimation->SetTargetReduction(DecimationTargetReduction);
            decimation->SetVolumePreservation(DecimationVolumePreservation);
            decimation->SetAttributeErrorMetric(DecimationAttributeErrorMetric);
            decimation->Update();
            if(errors->failed) throw std::runtime_error("ParaView output decimation failed");
            data = decimation->GetOutput();
        }
    }
    return data;
}
} // namespace
#endif
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
    using Clock=std::chrono::steady_clock;
    const auto begin=Clock::now();
    if(options.profile) *options.profile={};
    volume.validate();
    for(int n:volume.grid.cells) if(n<3) throw std::invalid_argument("ParaView AMR requires at least 3 cells per axis");
    if(!std::isfinite(options.isovalue)) throw std::invalid_argument("Invalid ParaView isovalue");
    for(double rho:volume.values) if(rho<0 || rho>1) throw std::invalid_argument("ParaView volume fractions must be in [0,1]");
    if(options.isovalue<0 || options.isovalue>1) throw std::invalid_argument("Volume fraction iso must be in [0,1]");
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
    const auto filterStart=Clock::now();
    contour->Update();
    const auto conversionStart=Clock::now();
    if(errors->failed) throw std::runtime_error("vtkAMRDualContour reported an error");
    auto* output=vtkMultiBlockDataSet::SafeDownCast(contour->GetOutputDataObject(0));
    if(!output) throw std::runtime_error("Unexpected ParaView output type");
    Mesh mesh;
    auto iterator=vtkSmartPointer<vtkCompositeDataIterator>::Take(output->NewIterator());
    iterator->SkipEmptyNodesOn();
    for(iterator->InitTraversal();!iterator->IsDoneWithTraversal();iterator->GoToNextItem()) {
        auto* poly=vtkPolyData::SafeDownCast(iterator->GetCurrentDataObject());
        if(!poly) continue;
        auto data=postprocessPolyData(poly,errors);
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
    if(options.profile) {
        options.profile->preparationMs=std::chrono::duration<double,std::milli>(filterStart-begin).count();
        options.profile->filterMs=std::chrono::duration<double,std::milli>(conversionStart-filterStart).count();
        // Includes optional smoothing/decimation, triangulation and Mesh conversion.
        options.profile->conversionMs=std::chrono::duration<double,std::milli>(Clock::now()-conversionStart).count();
    }
    return mesh;
#endif
}
}

