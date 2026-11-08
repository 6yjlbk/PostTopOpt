#include "DualContouring.h"
#include <Eigen/SVD>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace geometry {
namespace {
std::size_t product(const std::array<int,3>& n) {
    std::size_t result = 1;
    for (int v : n) {
        if (v <= 0 || result > std::numeric_limits<std::size_t>::max()/static_cast<std::size_t>(v))
            throw std::invalid_argument("Invalid or overflowing grid dimensions");
        result *= static_cast<std::size_t>(v);
    }
    return result;
}
// Solve relative to the intersection centroid to stabilize rank-deficient cells.
Vec3 place(const std::vector<Vec3>& points, const std::vector<Vec3>& normals,
           const Vec3& lo, const Vec3& hi, const ContourOptions& options) {
    Vec3 center = Vec3::Zero();
    for (const auto& p : points) center += p;
    center /= static_cast<double>(points.size());
    const auto count = static_cast<Eigen::Index>(points.size());
    Eigen::MatrixXd a(count + 3, 3);
    Eigen::VectorXd b(count + 3);
    for (Eigen::Index i=0; i<count; ++i) {
        a.row(i) = normals[i].transpose();
        b(i) = normals[i].dot(points[i]-center);
    }
    a.bottomRows(3) = std::sqrt(options.regularization)*Eigen::Matrix3d::Identity();
    b.tail(3).setZero();
    Vec3 best = center;
    double bestError = std::numeric_limits<double>::infinity();
    // Each axis is free, fixed at its lower bound, or fixed at its upper bound.
    // Enumerating 27 active sets covers the interior, faces, edges and corners.
    for (int state=0; state<(options.constrainToCell ? 27 : 1); ++state) {
        int code=state;
        std::vector<int> freeAxes;
        Vec3 offset=Vec3::Zero();
        for (int d=0; d<3; ++d) {
            int mode=code%3; code/=3;
            if (mode==0) freeAxes.push_back(d);
            else offset[d]=(mode==1 ? lo[d] : hi[d])-center[d];
        }
        if (!freeAxes.empty()) {
            Eigen::MatrixXd reduced(a.rows(), freeAxes.size());
            for (std::size_t j=0; j<freeAxes.size(); ++j) reduced.col(j)=a.col(freeAxes[j]);
            Eigen::VectorXd rhs=b-a*offset;
            Eigen::VectorXd sol=reduced.jacobiSvd(Eigen::ComputeThinU|Eigen::ComputeThinV).solve(rhs);
            for (std::size_t j=0; j<freeAxes.size(); ++j) offset[freeAxes[j]]=sol[j];
        }
        Vec3 candidate=center+offset;
        if (!candidate.allFinite()) continue;
        if (options.constrainToCell && ((candidate.array()<lo.array()-1e-12).any() ||
                                        (candidate.array()>hi.array()+1e-12).any())) continue;
        double error=(a*offset-b).squaredNorm();
        if (error<bestError) { bestError=error; best=candidate; }
    }
    if (!std::isfinite(bestError)) throw std::runtime_error("QEF failed");
    return options.constrainToCell ? best.cwiseMax(lo).cwiseMin(hi) : best;
}
}
Mesh dualContour(const Grid& grid, const Field& field, const Gradient& gradient,
                 const ContourOptions& options) {
    if (!field || !grid.origin.allFinite() || !grid.spacing.allFinite() ||
        (grid.spacing.array()<=0).any() || !std::isfinite(options.isovalue) ||
        !std::isfinite(options.regularization) || options.regularization<0 ||
        !std::isfinite(options.gradientStep) || options.gradientStep<=0)
        throw std::invalid_argument("Invalid field, grid or options");
    product(grid.cells);
    auto nodes=grid.cells;
    for (auto& n : nodes) {
        if (n==std::numeric_limits<int>::max()) throw std::invalid_argument("Grid too large");
        ++n;
    }
    const auto index=[](int x,int y,int z,const std::array<int,3>& n) {
        return (static_cast<std::size_t>(x)*n[1]+y)*n[2]+z;
    };
    const auto position=[&](int x,int y,int z)->Vec3 {
        return grid.origin+grid.spacing.cwiseProduct(Vec3(x,y,z));
    };
    const auto evaluate=[&](const Vec3& p) {
        double value=field(p)-options.isovalue;
        if (!std::isfinite(value)) throw std::invalid_argument("Non-finite field value");
        return value;
    };
    std::vector<double> values(product(nodes));
    for(int x=0;x<nodes[0];++x) for(int y=0;y<nodes[1];++y) for(int z=0;z<nodes[2];++z)
        values[index(x,y,z,nodes)]=evaluate(position(x,y,z));
    const auto value=[&](int x,int y,int z) { return values[index(x,y,z,nodes)]; };
    const auto missing=std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> ids(product(grid.cells),missing);
    Mesh mesh;
    for(int x=0;x<grid.cells[0];++x) for(int y=0;y<grid.cells[1];++y) for(int z=0;z<grid.cells[2];++z) {
        std::vector<Vec3> points,normals;
        for(int axis=0;axis<3;++axis) for(int u=0;u<2;++u) for(int v=0;v<2;++v) {
            std::array<int,3> start{x,y,z};
            start[(axis+1)%3]+=u; start[(axis+2)%3]+=v;
            auto end=start; ++end[axis];
            double f0=value(start[0],start[1],start[2]), f1=value(end[0],end[1],end[2]);
            if ((f0<0)==(f1<0)) continue;
            // Scaling avoids overflow when opposite finite values have a huge magnitude.
            double scale=std::max(std::abs(f0),std::abs(f1));
            double t=(f0/scale)/((f0/scale)-(f1/scale));
            Vec3 p=position(start[0],start[1],start[2]); p[axis]+=t*grid.spacing[axis];
            Vec3 n;
            if (gradient) n=gradient(p);
            else for(int d=0;d<3;++d) {
                Vec3 delta=Vec3::Zero(); delta[d]=options.gradientStep;
                n[d]=(evaluate(p+delta)-evaluate(p-delta))/(2*options.gradientStep);
            }
            const double length=n.stableNorm();
            if (!n.allFinite() || !std::isfinite(length) || length<=1e-14)
                throw std::invalid_argument("Undefined normal at an intersection");
            points.push_back(p); normals.push_back(n/length);
        }
        if(points.empty()) continue;
        ids[index(x,y,z,grid.cells)]=mesh.vertices.size();
        Vec3 lo=position(x,y,z);
        mesh.vertices.push_back(place(points,normals,lo,lo+grid.spacing,options));
    }
    // Cyclic transverse axes give a +axis orientation. Reverse when outside is at the start.
    for(int axis=0;axis<3;++axis) {
        int u=(axis+1)%3,v=(axis+2)%3;
        for(int k=0;k<grid.cells[axis];++k) for(int i=1;i<grid.cells[u];++i) for(int j=1;j<grid.cells[v];++j) {
            std::array<int,3> start{}; start[axis]=k; start[u]=i; start[v]=j;
            auto end=start; ++end[axis];
            double a=value(start[0],start[1],start[2]),b=value(end[0],end[1],end[2]);
            if((a<0)==(b<0)) continue;
            std::array<std::size_t,4> q{};
            constexpr int du[4]={-1,0,0,-1},dv[4]={-1,-1,0,0};
            for(int c=0;c<4;++c) {
                auto cell=start; cell[u]+=du[c]; cell[v]+=dv[c];
                q[c]=ids[index(cell[0],cell[1],cell[2],grid.cells)];
                if(q[c]==missing) throw std::runtime_error("Missing active cell");
            }
            if(a>=0) std::reverse(q.begin(),q.end());
            mesh.quads.push_back(q);
        }
    }
    return mesh;
}
}
