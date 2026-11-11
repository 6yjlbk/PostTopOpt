#include "CgalExtraction.h"
#include <cmath>
#include <type_traits>
#include <limits>
#include <stdexcept>
#ifdef POSTTOPOPT_WITH_CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/Octree_partition.h>
#include <CGAL/Isosurfacing_3/Value_function_3.h>
#include <CGAL/Isosurfacing_3/Gradient_function_3.h>
#include <CGAL/Isosurfacing_3/Dual_contouring_domain_3.h>
#include <CGAL/Isosurfacing_3/dual_contouring_3.h>
#endif
namespace geometry {
bool cgalAvailable() noexcept {
#ifdef POSTTOPOPT_WITH_CGAL
    return true;
#else
    return false;
#endif
}
Mesh cgalDualContour(const Grid& grid, const Field& field, const Gradient& gradient,
                     const CgalContourOptions& options) {
#ifndef POSTTOPOPT_WITH_CGAL
    (void)grid; (void)field; (void)gradient; (void)options;
    throw std::runtime_error("CGAL is disabled. Configure with -DPOSTTOPOPT_WITH_CGAL=ON and rebuild.");
#else
    if (!field || !grid.origin.allFinite() || !grid.spacing.allFinite() ||
        (grid.spacing.array()<=0).any() || !std::isfinite(options.isovalue) ||
        !std::isfinite(options.gradientStep) || options.gradientStep<=0)
        throw std::invalid_argument("Invalid CGAL field, grid or options");
    std::array<std::size_t,3> dimensions{};
    std::size_t count=1;
    for (int d=0;d<3;++d) {
        if (grid.cells[d]<=0 || grid.cells[d]==std::numeric_limits<int>::max())
            throw std::invalid_argument("Invalid CGAL grid dimensions");
        dimensions[d]=static_cast<std::size_t>(grid.cells[d])+1;
        if (count>std::numeric_limits<std::size_t>::max()/dimensions[d])
            throw std::invalid_argument("CGAL grid size overflow");
        count*=dimensions[d];
    }
    const Vec3 hi=grid.origin+grid.spacing.cwiseProduct(
        Vec3(grid.cells[0],grid.cells[1],grid.cells[2]));
    if (!hi.allFinite() || (hi.array()<=grid.origin.array()).any())
        throw std::invalid_argument("Invalid CGAL grid bounds");
    using Kernel=CGAL::Simple_cartesian<double>;
    using Point=Kernel::Point_3;
    using Vector=Kernel::Vector_3;
    using CGrid=CGAL::Isosurfacing::Cartesian_grid_3<Kernel>;
    const auto evaluate=[&](const Vec3& p) {
        const double v=field(p);
        if (!std::isfinite(v) || !std::isfinite(v-options.isovalue))
            throw std::invalid_argument("Non-finite CGAL field value");
        return v;
    };
    const auto value=[&](const Point& p) {return evaluate(Vec3(p.x(),p.y(),p.z()));};
    const auto normal=[&](const Point& point)->Vector {
        const Vec3 p(point.x(),point.y(),point.z());
        Vec3 g;
        if (gradient) g=gradient(p);
        else for (int d=0;d<3;++d) {
            Vec3 delta=Vec3::Zero(); delta[d]=options.gradientStep;
            g[d]=(evaluate(p+delta)-evaluate(p-delta))/(2*options.gradientStep);
        }
        const double length=g.stableNorm();
        if (!g.allFinite() || !std::isfinite(length) || length<=1e-14)
            throw std::invalid_argument("Undefined CGAL normal");
        g/=length;
        return Vector(g.x(),g.y(),g.z());
    };
    std::vector<Point> points;
    std::vector<std::vector<std::size_t>> polygons;
    const auto extractPartition = [&](const auto& partition) {
        using Partition=std::decay_t<decltype(partition)>;
        CGAL::Isosurfacing::Value_function_3<Partition> values(value,partition);
        CGAL::Isosurfacing::Gradient_function_3<Partition> gradients(normal,partition);
        const auto extract = [&](const auto& intersectionOracle) {
            const auto domain=CGAL::Isosurfacing::create_dual_contouring_domain_3(
                partition,values,gradients,intersectionOracle);
            CGAL::Isosurfacing::dual_contouring<CGAL::Sequential_tag>(
                domain,options.isovalue,points,polygons,
                CGAL::parameters::do_not_triangulate_faces(true).constrain_to_cell(options.constrainToCell));
        };
        if (options.bisection)
            extract(CGAL::Isosurfacing::Dichotomy_edge_intersection{}); // Биссекция.
        else
            extract(CGAL::Isosurfacing::Linear_interpolation_edge_intersection{});
    };
    if (options.octree) {
        if (options.octreeMinDepth>options.octreeMaxDepth || options.octreeMaxDepth>20)
            throw std::invalid_argument("Octree requires minDepth <= maxDepth <= 20");
        using Octree=CGAL::Octree<Kernel,std::vector<Point>>;
        std::vector<Point> bounds{
            Point(grid.origin.x(),grid.origin.y(),grid.origin.z()),
            Point(hi.x(),hi.y(),hi.z())};
        Octree partition(bounds);
        // Experimental sign-based refinement, as in CGAL's contouring_octree example.
        // Equal-sign corners can hide a surface; mixed depths can produce cracks.
        partition.refine([&](const Octree::Node_index& node,const Octree& tree) {
            if (tree.depth(node)<options.octreeMinDepth) return true;
            if (tree.depth(node)>=options.octreeMaxDepth) return false;
            const auto box=tree.bbox(node);
            bool inside=false,outside=false;
            for (double x : {box.xmin(),box.xmax()})
                for (double y : {box.ymin(),box.ymax()})
                    for (double z : {box.zmin(),box.zmax()}) {
                        const double v=value(Point(x,y,z));
                        inside |= v<=options.isovalue;
                        outside |= v>options.isovalue;
                    }
            return inside && outside;
        });
        extractPartition(partition);
    } else {
        // CGAL expects NODE counts, whereas our Grid stores CELL counts.
        CGrid partition(Point(grid.origin.x(),grid.origin.y(),grid.origin.z()),
                        Point(hi.x(),hi.y(),hi.z()),dimensions);
        extractPartition(partition);
    }
    Mesh mesh;
    mesh.vertices.reserve(points.size());
    for (const auto& p:points) {
        Vec3 v(p.x(),p.y(),p.z());
        if (!v.allFinite()) throw std::runtime_error("CGAL returned a non-finite vertex");
        mesh.vertices.push_back(v);
    }
    for (const auto& face:polygons) {
        for (auto i:face) if (i>=mesh.vertices.size())
            throw std::runtime_error("CGAL returned an invalid index");
        if (face.size()==4) mesh.quads.push_back({face[0],face[1],face[2],face[3]});
        else if (face.size()==3) mesh.triangles.push_back({face[0],face[1],face[2]});
        else throw std::runtime_error("Unexpected CGAL polygon size");
    }
    return mesh;
#endif
}
}
