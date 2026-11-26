#include "DualContouring.h"
#include <Eigen/LU>
#include <Eigen/Geometry>
#include <functional>
#include <utility>
#include <Eigen/SVD>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace postopdc {
void ElementField::shape(const Vec3& q, std::array<float,8>& n,
                         std::array<Vec3,8>& d) const {
    const float r=q.x(),s=q.y(),t=q.z();
    if (element->type==ElementType::Tetra4) {
        n={1-r-s-t,r,s,t,0,0,0,0};
        d[0]=Vec3(-1,-1,-1); d[1]=Vec3(1,0,0);
        d[2]=Vec3(0,1,0); d[3]=Vec3(0,0,1);
    } else if (element->type==ElementType::Wedge6) {
        const float b[3]={1-r-s,r,s};
        const Vec3 db[3]={Vec3(-1,-1,0),Vec3(1,0,0),Vec3(0,1,0)};
        for (int a=0;a<3;++a) {
            n[a]=b[a]*(1-t); n[a+3]=b[a]*t;
            d[a]=db[a]*(1-t); d[a].z()=-b[a];
            d[a+3]=db[a]*t; d[a+3].z()=b[a];
        }
    } else {
        static constexpr int corner[8][3]={
            {0,0,0},{1,0,0},{1,1,0},{0,1,0},
            {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
        for (int a=0;a<8;++a) {
            const float x=corner[a][0]?r:1-r, y=corner[a][1]?s:1-s,
                        z=corner[a][2]?t:1-t;
            n[a]=x*y*z;
            d[a]=Vec3((corner[a][0]?1.f:-1.f)*y*z,
                      (corner[a][1]?1.f:-1.f)*x*z,
                      (corner[a][2]?1.f:-1.f)*x*y);
        }
    }
}
Vec3 ElementField::position(const Vec3& q) const {
    std::array<float,8> n; std::array<Vec3,8> d; shape(q,n,d);
    // Relative coordinates reduce cancellation at large world offsets.
    const Vec3 origin=points[element->nodes[0]];
    Vec3 p=origin;
    for(int a=1;a<element->nodeCount();++a)
        p+=n[a]*(points[element->nodes[a]]-origin);
    return p;
}
float ElementField::value(const Vec3& q) const {
    std::array<float,8> n; std::array<Vec3,8> d; shape(q,n,d);
    float v=0;
    for(int a=0;a<element->nodeCount();++a) v+=n[a]*nodalVF[element->nodes[a]];
    return v;
}
Eigen::Matrix3f ElementField::jacobian(const Vec3& q) const {
    std::array<float,8> n; std::array<Vec3,8> d; shape(q,n,d);
    Eigen::Matrix3f j=Eigen::Matrix3f::Zero();
    const Vec3 origin=points[element->nodes[0]];
    for(int a=1;a<element->nodeCount();++a)
        j+=(points[element->nodes[a]]-origin)*d[a].transpose();
    return j;
}
Vec3 ElementField::gradient(const Vec3& q) const {
    std::array<float,8> n; std::array<Vec3,8> d; shape(q,n,d);
    Vec3 g=Vec3::Zero();
    for(int a=1;a<element->nodeCount();++a)
        g+=d[a]*(nodalVF[element->nodes[a]]-nodalVF[element->nodes[0]]);
    const Eigen::Matrix3f jt=jacobian(q).transpose();
    const auto lu=jt.fullPivLu();
    if(!lu.isInvertible()) throw std::runtime_error("Singular FE Jacobian");
    return lu.solve(g);
}
Vec3 ElementField::project(const Vec3& q) const {
    Vec3 p=q.cwiseMax(0.f);
    if(element->type==ElementType::Hexa8) return p.cwiseMin(1.f);
    const int dimensions=element->type==ElementType::Tetra4?3:2;
    if(dimensions==2) p.z()=(std::min)(p.z(),1.f);
    if(p.head(dimensions).sum()<=1) return p;
    std::array<float,3> sorted{q.x(),q.y(),q.z()};
    std::sort(sorted.begin(),sorted.begin()+dimensions,std::greater<float>());
    float sum=0,theta=0;
    for(int i=0;i<dimensions;++i) {
        sum+=sorted[i];
        const float candidate=(sum-1)/(i+1);
        if(sorted[i]>candidate) theta=candidate;
    }
    for(int i=0;i<dimensions;++i) p[i]=(std::max)(0.f,q[i]-theta);
    return p;
}
namespace {
using DVec=Eigen::Vector3d;
using Key=std::array<std::size_t,4>;
constexpr auto missing=(std::numeric_limits<std::size_t>::max)();
struct Hash {
    std::size_t operator()(const Key& k) const {
        std::size_t h=0;
        for(auto v:k) h^=v+0x9e3779b9+(h<<6)+(h>>2);
        return h;
    }
};
Key edgeKey(std::size_t a,std::size_t b) {
    if(a>b) std::swap(a,b);
    return {a,b,missing,missing};
}
struct Topology {
    std::vector<Vec3> corners;
    std::vector<std::array<int,2>> edges;
    std::vector<std::vector<int>> faces;
};
const Topology& topology(ElementType type) {
    static const Topology tetra{
        {{0,0,0},{1,0,0},{0,1,0},{0,0,1}},
        {{0,1},{1,2},{2,0},{0,3},{1,3},{2,3}},
        {{0,2,1},{0,1,3},{1,2,3},{2,0,3}}};
    static const Topology wedge{
        {{0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,0,1},{0,1,1}},
        {{0,1},{1,2},{2,0},{3,4},{4,5},{5,3},{0,3},{1,4},{2,5}},
        {{0,2,1},{3,4,5},{0,1,4,3},{1,2,5,4},{2,0,3,5}}};
    static const Topology hexa{
        {{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}},
        {{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}},
        {{0,3,2,1},{4,5,6,7},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7}}};
    return type==ElementType::Tetra4?tetra:type==ElementType::Wedge6?wedge:hexa;
}
// Double precision QEF, centered at the mean Hermite intersection.
// On mapped elements the constrained problem is nonlinear. Projected descent
// gives a feasible local solution, not a certified global minimum.
Vec3 place(const ElementField& field, const std::array<Vec3,12>& p,
           const std::array<Vec3,12>& n,int count,Vec3 q,
           const ExtractionOptions& options) {
    DVec center=DVec::Zero();
    for(int i=0;i<count;++i) center+=p[i].cast<double>();
    center/=count;
    Eigen::Matrix<double,Eigen::Dynamic,3> a(count+3,3);
    Eigen::VectorXd b(count+3);
    for(int i=0;i<count;++i) {
        a.row(i)=n[i].cast<double>().transpose();
        b[i]=n[i].cast<double>().dot(p[i].cast<double>()-center);
    }
    a.bottomRows(3)=std::sqrt(options.regularizationEnabled?
        double(options.regularization):0.)*Eigen::Matrix3d::Identity();
    b.tail(3).setZero();
    const DVec target=center+a.jacobiSvd(Eigen::ComputeFullU|Eigen::ComputeFullV).solve(b);
    if(!options.constrainToCell) return target.cast<float>();
    // Fast path: invert the FE map for the unconstrained QEF minimum.
    // Accept only an interior, converged inverse; never clamp a world-space AABB.
    Vec3 inverse=q;
    for(int iteration=0;iteration<20;++iteration) {
        const Eigen::Matrix3d j=field.jacobian(inverse).cast<double>();
        const auto lu=j.fullPivLu();
        if(!lu.isInvertible()) break;
        const DVec residual=target-field.position(inverse).cast<double>();
        if(residual.norm()<=1e-6*j.norm() &&
           (field.project(inverse)-inverse).squaredNorm()<=1e-14f)
            return field.position(field.project(inverse));
        const DVec delta=lu.solve(residual);
        if(!delta.allFinite()) break;
        inverse+=delta.cast<float>();
        if(!inverse.allFinite() || inverse.norm()>1e4f) break;
    }
    const auto error=[&](const Vec3& at) {
        return (a*(field.position(at).cast<double>()-center)-b).squaredNorm();
    };
    q=field.project(q);
    double current=error(q);
    for(int iteration=0;iteration<128;++iteration) {
        const DVec residualPosition=field.position(q).cast<double>()-center;
        const Eigen::Matrix3d j=field.jacobian(q).cast<double>();
        const Eigen::MatrixXd aj=a*j;
        const DVec g=aj.transpose()*(a*residualPosition-b);
        double step=1./(std::max)(aj.squaredNorm(),1e-30);
        bool accepted=false;
        for(int search=0;search<24;++search) {
            const Vec3 next=field.project((q.cast<double>()-step*g).cast<float>());
            if((next-q).squaredNorm()<1e-14f) return field.position(q);
            const double candidate=error(next);
            if(candidate<current) {q=next;current=candidate;accepted=true;break;}
            step*=0.5;
        }
        if(!accepted) break;
    }
    return field.position(q);
}
struct Face {
    std::vector<std::size_t> nodes;
    std::array<std::size_t,2> cells{missing,missing};
};
struct Edge {
    float t=0;
    bool boundary=false;
    std::vector<std::array<std::size_t,2>> links;
};
} // namespace
Surface dualContour(const std::vector<Vec3>& points,
    const std::vector<const Element*>& elements,const std::vector<float>& nodalVF,
    const ExtractionOptions& options) {
    Surface result;
    std::vector<std::size_t> vertex(elements.size(),missing);
    std::unordered_map<Key,Face,Hash> faces;
    std::unordered_map<Key,Edge,Hash> edges;
    std::vector<Key> edgeOrder; // Stable output independent of hash iteration order.
    for(std::size_t cell=0;cell<elements.size();++cell) {
        const Element* element=elements[cell];
        const auto& topo=topology(element->type);
        ElementField field(element,points,nodalVF);
        for(const auto& local:topo.faces) {
            Key key{missing,missing,missing,missing};
            for(std::size_t a=0;a<local.size();++a) key[a]=element->nodes[local[a]];
            std::sort(key.begin(),key.end());
            auto& face=faces[key];
            if(face.cells[0]==missing) {
                face.cells[0]=cell;
                for(int a:local) face.nodes.push_back(element->nodes[a]);
            } else {
                if(face.cells[1]!=missing) throw std::runtime_error("Nonmanifold FE face");
                face.cells[1]=cell;
            }
        }
        std::array<Vec3,12> samples,normals;
        Vec3 meanQ=Vec3::Zero();
        int count=0;
        for(const auto& local:topo.edges) {
            int i=local[0],j=local[1];
            if(element->nodes[i]>element->nodes[j]) std::swap(i,j);
            const auto key=edgeKey(element->nodes[i],element->nodes[j]);
            const float f0=options.iso-nodalVF[key[0]], f1=options.iso-nodalVF[key[1]];
            if((f0<0)==(f1<0)) continue;
            auto inserted=edges.try_emplace(key);
            Edge& edge=inserted.first->second;
            const Vec3 q0=topo.corners[i],q1=topo.corners[j];
            if(inserted.second) {
                edgeOrder.push_back(key);
                if(f0==0) edge.t=0;
                else if(f1==0) edge.t=1;
                else if(!options.bisection) edge.t=f0/(f0-f1);
                else {
                    float left=0,right=1;
                    for(int iteration=0;iteration<24;++iteration) {
                        edge.t=(left+right)*0.5f;
                        const float residual=options.iso-field.value(q0+edge.t*(q1-q0));
                        if(std::abs(residual)<=1e-7f) break;
                        if((residual<0)==(f0<0)) left=edge.t; else right=edge.t;
                    }
                }
            }
            const Vec3 q=q0+edge.t*(q1-q0);
            // Reuse the intersection only: gradients may jump across FE faces.
            const Vec3 gradient=-field.gradient(q);
            const float length=gradient.stableNorm();
            if(!std::isfinite(length) || length<=0)
                throw std::runtime_error("Undefined Hermite normal");
            samples[count]=field.position(q); normals[count]=gradient/length;
            meanQ+=q; ++count;
        }
        if(!count) continue;
        const Vec3 p=place(field,samples,normals,count,meanQ/float(count),options);
        vertex[cell]=result.vertices.size();
        result.vertices.push_back({p.x(),p.y(),p.z()});
    }
    // Shared faces determine the cyclic order around each intersected edge.
    for(const auto& item:faces) {
        const auto& face=item.second;
        for(std::size_t i=0;i<face.nodes.size();++i) {
            auto found=edges.find(edgeKey(face.nodes[i],face.nodes[(i+1)%face.nodes.size()]));
            if(found==edges.end()) continue;
            if(face.cells[1]==missing) found->second.boundary=true;
            else found->second.links.push_back(face.cells);
        }
    }
    for(const auto& key:edgeOrder) {
        const auto& edge=edges.at(key);
        if(edge.boundary) continue; // No invented exterior elements or caps.
        const auto& links=edge.links;
        if(links.size()<3) throw std::runtime_error("Invalid interior FE edge ring");
        std::size_t first=missing;
        for(const auto& link:links) first=(std::min)(first,(std::min)(link[0],link[1]));
        std::size_t previous=missing,current=first;
        std::vector<std::size_t> polygon;
        Vec3 area=Vec3::Zero(),origin=points[key[0]];
        std::vector<Vec3> centers;
        for(std::size_t step=0;step<links.size();++step) {
            if(vertex[current]==missing) throw std::runtime_error("Missing active FE vertex");
            polygon.push_back(vertex[current]);
            const Element* e=elements[current];
            Vec3 center=Vec3::Zero();
            for(int a=0;a<e->nodeCount();++a) center+=points[e->nodes[a]]-origin;
            centers.push_back(center/float(e->nodeCount()));
            std::size_t next=missing; int degree=0;
            for(const auto& link:links) {
                if(link[0]!=current && link[1]!=current) continue;
                ++degree;
                const auto other=link[0]==current?link[1]:link[0];
                if(other!=previous) next=(std::min)(next,other);
            }
            if(degree!=2 || next==missing) throw std::runtime_error("Nonmanifold FE edge ring");
            previous=current; current=next;
            if(current==first && step+1!=links.size())
                throw std::runtime_error("Disconnected FE edge ring");
        }
        if(current!=first) throw std::runtime_error("Open interior FE edge ring");
        for(std::size_t i=0;i<centers.size();++i)
            area+=centers[i].cross(centers[(i+1)%centers.size()]);
        Vec3 outward=points[key[1]]-points[key[0]];
        if(nodalVF[key[0]]<=options.iso) outward=-outward;
        if(area.dot(outward)<0) std::reverse(polygon.begin()+1,polygon.end());
        if(polygon.size()==3) result.triangles.push_back({polygon[0],polygon[1],polygon[2]});
        else if(polygon.size()==4) result.quads.push_back({polygon[0],polygon[1],polygon[2],polygon[3]});
        else result.polygons.push_back(std::move(polygon));
    }
    return result;
}
} // namespace postopdc

