#include "PlyDistanceField.h"
#include <vtkCellArray.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkSmartPointer.h>
#include <Eigen/Geometry>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace geometry {
PlyDistanceField makePlyDistanceField(const std::filesystem::path& path, int resolution, int padding) {
    if (resolution<4 || resolution>512 || padding<2 || padding>32)
        throw std::invalid_argument("Resolution must be 4..512; padding must be 2..32 cells");
    const auto imported=readPly(path);
    if (!imported.report.closedByEdgeIncidence())
        throw std::invalid_argument("Signed-distance extraction requires a closed input without non-manifold edges; use PLY import diagnostics first");
    vtkNew<vtkPoints> points;
    points->SetDataTypeToDouble();
    for (const auto& p:imported.mesh.vertices) points->InsertNextPoint(p.data());
    vtkNew<vtkCellArray> faces;
    for (const auto& face:imported.mesh.triangles) {
        const auto& a=imported.mesh.vertices[face[0]];
        const auto& b=imported.mesh.vertices[face[1]];
        const auto& c=imported.mesh.vertices[face[2]];
        if ((b-a).cross(c-a).squaredNorm()==0)
            throw std::invalid_argument("Degenerate source triangle: repair the input explicitly");
        vtkIdType ids[3]={static_cast<vtkIdType>(face[0]),static_cast<vtkIdType>(face[1]),static_cast<vtkIdType>(face[2])};
        faces->InsertNextCell(3,ids);
    }
    vtkNew<vtkPolyData> poly;
    poly->SetPoints(points); poly->SetPolys(faces);
    // Prepare orientation on an in-memory copy; do not alter the source PLY.
    vtkNew<vtkPolyDataNormals> normals;
    normals->SetInputData(poly);
    normals->ConsistencyOn(); normals->AutoOrientNormalsOn();
    normals->SplittingOff(); normals->NonManifoldTraversalOff();
    normals->ComputePointNormalsOff(); normals->ComputeCellNormalsOn();
    normals->Update();
    if (!normals->GetOutput() || normals->GetOutput()->GetNumberOfPolys()!=faces->GetNumberOfCells())
        throw std::runtime_error("Surface orientation failed");
    auto distance=vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
    distance->SetInput(normals->GetOutput());
    const Vec3 extent=imported.report.maximum-imported.report.minimum;
    const double h=extent.maxCoeff()/resolution;
    if (!std::isfinite(h) || h<=0) throw std::invalid_argument("Invalid model bounds");
    PlyDistanceField result;
    result.source=imported.report;
    result.grid.origin=imported.report.minimum-Vec3::Constant(padding*h);
    result.grid.spacing=Vec3::Constant(h);
    for (int d=0;d<3;++d) result.grid.cells[d]=static_cast<int>(std::ceil(extent[d]/h))+2*padding;
    result.field=[distance](const Vec3& p) {
        double xyz[3]={p.x(),p.y(),p.z()};
        double value=distance->EvaluateFunction(xyz);
        if (!std::isfinite(value)) throw std::runtime_error("Non-finite VTK distance");
        return value;
    };
    // Detect obvious inverted/outside signs; not a complete topology validation.
    for (int corner=0;corner<8;++corner) {
        Vec3 p=result.grid.origin;
        for(int d=0;d<3;++d) if(corner&(1<<d)) p[d]+=result.grid.spacing[d]*result.grid.cells[d];
        if(result.field(p)<=0) throw std::runtime_error("Distance sign outside model is not positive");
    }
    return result;
}
}
