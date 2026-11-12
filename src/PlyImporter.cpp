#include "PlyImporter.h"
#include <vtkCellArray.h>
#include <vtkCommand.h>
#include <vtkFeatureEdges.h>
#include <vtkNew.h>
#include <vtkPLYReader.h>
#include <vtkPolyData.h>
#include <vtkTriangleFilter.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace geometry {
namespace {
class ErrorObserver : public vtkCommand {
public:
    static ErrorObserver* New() { return new ErrorObserver; }
    bool failed = false;
    std::string message;
    void Execute(vtkObject*, unsigned long, void* data) override {
        failed = true;
        message = data ? static_cast<const char*>(data) : "Unknown VTK error";
    }
};
std::size_t countEdges(vtkPolyData* data, bool boundary) {
    vtkNew<vtkFeatureEdges> filter;
    vtkNew<ErrorObserver> errors;
    filter->AddObserver(vtkCommand::ErrorEvent, errors);
    filter->SetInputData(data);
    filter->BoundaryEdgesOff();
    filter->NonManifoldEdgesOff();
    filter->FeatureEdgesOff();
    filter->ManifoldEdgesOff();
    if (boundary) filter->BoundaryEdgesOn();
    else filter->NonManifoldEdgesOn();
    filter->Update();
    if (errors->failed) throw std::runtime_error(errors->message);
    return static_cast<std::size_t>(filter->GetOutput()->GetNumberOfLines());
}
}
ImportedPly readPly(const std::filesystem::path& path) {
    if (!std::filesystem::is_regular_file(path))
        throw std::invalid_argument("PLY file does not exist: " + path.string());
    vtkNew<ErrorObserver> errors;
    vtkNew<vtkPLYReader> reader;
    reader->AddObserver(vtkCommand::ErrorEvent, errors);
    const auto filename = path.string();
    reader->SetFileName(filename.c_str());
    reader->SetDuplicatePointsForFaceTexture(false);
    reader->Update();
    if (errors->failed) throw std::runtime_error("PLY reader: " + errors->message);
    auto* source = reader->GetOutput();
    if (!source || source->GetNumberOfPoints() == 0 || source->GetNumberOfPolys() == 0)
        throw std::runtime_error("PLY must contain vertices and polygon faces");
    if (source->GetNumberOfVerts() || source->GetNumberOfLines() || source->GetNumberOfStrips())
        throw std::runtime_error("PLY contains unsupported non-polygon cells");
    ImportedPly result;
    result.report.sourcePoints = static_cast<std::size_t>(source->GetNumberOfPoints());
    result.report.sourcePolygons = static_cast<std::size_t>(source->GetNumberOfPolys());
    vtkNew<vtkTriangleFilter> triangles;
    triangles->AddObserver(vtkCommand::ErrorEvent, errors);
    triangles->SetInputData(source);
    triangles->PassVertsOff();
    triangles->PassLinesOff();
    triangles->Update();
    if (errors->failed) throw std::runtime_error("PLY triangulation: " + errors->message);
    auto* data = triangles->GetOutput();
    if (!data || data->GetNumberOfPolys() == 0) throw std::runtime_error("Empty triangulated PLY");
    result.report.minimum = Vec3::Constant(std::numeric_limits<double>::infinity());
    result.report.maximum = -result.report.minimum;
    result.mesh.vertices.reserve(static_cast<std::size_t>(data->GetNumberOfPoints()));
    for (vtkIdType i=0; i<data->GetNumberOfPoints(); ++i) {
        double xyz[3]; data->GetPoint(i, xyz);
        Vec3 point(xyz[0],xyz[1],xyz[2]);
        if (!point.allFinite()) throw std::runtime_error("Non-finite PLY coordinates");
        result.report.minimum = result.report.minimum.cwiseMin(point);
        result.report.maximum = result.report.maximum.cwiseMax(point);
        result.mesh.vertices.push_back(point);
    }
    vtkIdType count = 0;
    const vtkIdType* ids = nullptr;
    auto* cells = data->GetPolys();
    cells->InitTraversal();
    while (cells->GetNextCell(count, ids)) {
        if (count != 3) throw std::runtime_error("Triangulation returned a non-triangle");
        std::array<std::size_t,3> face{};
        for (int j=0;j<3;++j) {
            if (ids[j]<0 || ids[j]>=data->GetNumberOfPoints()) throw std::runtime_error("Invalid PLY index");
            face[j]=static_cast<std::size_t>(ids[j]);
        }
        if (face[0]==face[1] || face[1]==face[2] || face[2]==face[0])
            throw std::runtime_error("PLY contains a triangle with repeated indices");
        result.mesh.triangles.push_back(face);
    }
    result.report.boundaryEdges = countEdges(data,true);
    result.report.nonManifoldEdges = countEdges(data,false);
    return result;
}
}
