#include "SurfaceExtraction.h"
#include "ObjWriter.h"
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
using namespace geometry;
void check(bool condition,const char* message) {
    if (!condition) throw std::runtime_error(message);
}
template<class F> void rejects(F action) {
    bool rejected=false;
    try { action(); } catch (const std::invalid_argument&) { rejected=true; }
    check(rejected,"Invalid input was accepted");
}
void checkMesh(const Mesh& mesh) {
    for (const auto& p:mesh.vertices) check(p.allFinite(),"Non-finite vertex");
    for (const auto& t:mesh.triangles) {
        for(auto id:t) check(id<mesh.vertices.size(),"Invalid triangle index");
        check(t[0]!=t[1] && t[1]!=t[2] && t[0]!=t[2],"Repeated triangle index");
        check((mesh.vertices[t[1]]-mesh.vertices[t[0]]).cross(
            mesh.vertices[t[2]]-mesh.vertices[t[0]]).squaredNorm()>0,"Collapsed triangle");
    }
}
int main() {
    try {
        // MC coverage: all eight-bit sign cases, interpolated within one cube.
        constexpr int c[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                               {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
        constexpr int e[12][2]={{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},
                                {6,7},{7,4},{0,4},{1,5},{2,6},{3,7}};
        Grid single; single.origin=Vec3::Zero(); single.spacing=Vec3::Ones(); single.cells={1,1,1};
        for (int mask=0;mask<256;++mask) {
            const auto f=[&](const Vec3& p) {
                double value=0;
                for(int i=0;i<8;++i) {
                    double weight=1;
                    for(int d=0;d<3;++d) weight*=c[i][d] ? p[d] : 1-p[d];
                    value+=((mask&(1<<i)) ? -1. : 1.)*weight;
                }
                return value;
            };
            const auto mesh=marchingCubes(single,f);
            checkMesh(mesh);
            check(mesh.quads.empty(),"MC produced quads");
            std::size_t crossings=0;
            for(const auto& edge:e) if (((mask>>edge[0])&1)!=((mask>>edge[1])&1)) ++crossings;
            check(mesh.vertices.size()==crossings,"MC case lost or duplicated an edge intersection");
            check(mesh.triangles.size()<=5,"Classic MC case exceeded five triangles");
            for(const auto& p:mesh.vertices) check(std::abs(f(p))<1e-12,"MC vertex is off an edge root");
        }
        // Non-cubic grid, nonzero isovalue, sharing across cells, outward winding.
        Grid g; g.origin=Vec3::Zero(); g.spacing=Vec3(0.2,0.3,0.4); g.cells={4,3,2};
        for(double iso:{0.35,0.4}) {
            const auto mesh=marchingCubes(g,[](const Vec3& p){return p.x();},iso);
            checkMesh(mesh);
            check(mesh.vertices.size()==12 && mesh.triangles.size()==12,"Plane not welded across cells");
            std::map<std::pair<std::size_t,std::size_t>,int> edgeUse;
            for(const auto& p:mesh.vertices) check(std::abs(p.x()-iso)<1e-12,"Wrong plane isovalue");
            for(const auto& t:mesh.triangles) {
                check((mesh.vertices[t[1]]-mesh.vertices[t[0]]).cross(
                    mesh.vertices[t[2]]-mesh.vertices[t[0]]).x()>0,"MC inward winding");
                for(int i=0;i<3;++i) ++edgeUse[std::minmax(t[i],t[(i+1)%3])];
            }
            int boundary=0;
            for(const auto& [edge,count]:edgeUse) {
                check(count==1 || count==2,"Non-manifold plane edge");
                boundary+=count==1;
            }
            check(boundary==10,"Unexpected interior cracks");
        }
        // Different data, same solver: an exact expected endpoint-plane result.
        Grid q; q.origin=Vec3(0.2,0,0); q.spacing=Vec3(0.6,1,1); q.cells={1,1,1};
        ContourOptions options; options.isovalue=0.25; options.regularization=0;
        const auto f=[](const Vec3& p){return p.x()*p.x();};
        const auto grad=[](const Vec3& p)->Vec3 {return Vec3(2*p.x(),0,0);};
        const auto dc=dualContour(q,f,grad,options);
        const auto emilk=dualContourEmilk(q,f,grad,options);
        check(dc.vertices.size()==1 && emilk.vertices.size()==1,"Missing DC vertex");
        check(std::abs(dc.vertices[0].x()-0.41)<1e-10,"Classic DC regression");
        check(std::abs(emilk.vertices[0].x()-0.640625)<1e-10,"Wrong endpoint-plane constraints");
        ContourOptions scaledOptions=options; scaledOptions.isovalue*=7;
        const auto scaled=dualContourEmilk(q,[](const Vec3& p){return 7*p.x()*p.x();},
            [](const Vec3& p)->Vec3 {return Vec3(14*p.x(),0,0);},scaledOptions);
        check((emilk.vertices[0]-scaled.vertices[0]).norm()<1e-10,"Field scaling changed normalized constraints");
        const auto numerical=dualContourEmilk(q,f,{},options);
        check((emilk.vertices[0]-numerical.vertices[0]).norm()<1e-8,"Numerical gradient disagrees");
        // Rank-deficient plane, same connectivity in both DC variants.
        options.isovalue=0.35;
        const auto plane=[](const Vec3& p){return p.x();};
        const auto a=dualContour(g,plane,{},options), b=dualContourEmilk(g,plane,{},options);
        check(a.quads==b.quads && !b.quads.empty(),"DC connectivity changed");
        for(const auto& p:b.vertices) check(std::abs(p.x()-0.35)<1e-10,"Emilk plane QEF failed");
        check(marchingCubes(g,[](const Vec3&){return 1.;}).vertices.empty(),"MC constant field");
        check(dualContourEmilk(g,[](const Vec3&){return 1.;}).vertices.empty(),"Emilk constant field");
        rejects([&]{Grid bad=g;bad.spacing.x()=0;marchingCubes(bad,plane);});
        rejects([&]{Grid bad=g;bad.cells[0]=-1;marchingCubes(bad,plane);});
        rejects([&]{marchingCubes(g,[](const Vec3&){return std::numeric_limits<double>::quiet_NaN();});});
        rejects([&]{dualContourEmilk(g,plane,[](const Vec3&)->Vec3{return Vec3::Zero();},options);});
        rejects([]{parseMethod("invalid");});
        // MC must never call the gradient.
        auto mc=extractSurface(ExtractionMethod::MarchingCubes,g,plane,
            [](const Vec3&)->Vec3 {throw std::runtime_error("MC queried gradient");},options);
        Mesh mixed=mc;
        mixed.quads.push_back({0,1,2,3});
        const auto path=std::filesystem::path("test-mixed.obj");
        writeObj(mixed,path);
        std::ifstream in(path); std::string line; std::size_t faces=0,triangles=0,quads=0;
        while(std::getline(in,line)) if(line.starts_with("f ")) {
            std::istringstream row(line.substr(2));std::size_t id,n=0;
            while(row>>id){check(id>0 && id<=mixed.vertices.size(),"OBJ index");++n;}
            ++faces;triangles+=n==3;quads+=n==4;
        }
        in.close();std::filesystem::remove(path);
        check(faces==mixed.triangles.size()+1 && triangles==mixed.triangles.size() && quads==1,"Mixed OBJ export failed");
        std::cout<<"256 MC cases, winding/sharing, exact-isovalue plane, corner QEF, scaling, validation and OBJ passed\n";
    } catch(const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
}

