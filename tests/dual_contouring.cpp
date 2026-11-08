#include <Eigen/Geometry>
#include <algorithm>
#include "DualContouring.h"
#include "ObjWriter.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
using namespace geometry;
void check(bool ok,const char* message) { if(!ok) throw std::runtime_error(message); }
int main() {
 try {
    Grid g; g.origin=Vec3(-1.7,-1.5,-1.3); g.spacing=Vec3(0.1,0.12,0.14); g.cells={34,25,19};
    auto sphere=[](const Vec3& p){return p.norm()-0.93;};
    auto mesh=dualContour(g,sphere);
    check(!mesh.quads.empty(),"Empty sphere");
    for(const auto& p:mesh.vertices) check(p.allFinite() && std::abs(sphere(p))<0.03,"Sphere approximation error");
    std::map<std::pair<std::size_t,std::size_t>,std::pair<int,int>> edges;
    for(const auto& q:mesh.quads) {
        Vec3 n=(mesh.vertices[q[1]]-mesh.vertices[q[0]]).cross(mesh.vertices[q[2]]-mesh.vertices[q[0]]);
        check(n.dot(mesh.vertices[q[0]])>0,"Incorrect face winding");
        for(int i=0;i<4;++i) {
            auto a=q[i],b=q[(i+1)%4]; check(a<mesh.vertices.size() && a!=b,"Invalid index");
            auto& e=edges[std::minmax(a,b)]; ++e.first; e.second+=(a<b ? 1 : -1);
        }
    }
    for(const auto& [key,e]:edges) check(e.first==2 && e.second==0,"Sphere is not closed and consistently oriented");
    check(static_cast<long long>(mesh.vertices.size())-static_cast<long long>(edges.size())+
          static_cast<long long>(mesh.quads.size())==2,"Unexpected sphere Euler characteristic");
    check(dualContour(g,[](const Vec3&){return 1.;}).vertices.empty(),"Constant field not empty");
    ContourOptions options; options.isovalue=0.25; options.regularization=0;
    auto plane=dualContour(g,[](const Vec3& p){return p.x();}, {},options);
    check(!plane.quads.empty(),"Empty plane");
    for(const auto& p:plane.vertices) check(std::abs(p.x()-0.25)<1e-10,"Rank deficient QEF failed");
    bool rejected=false;
    try { Grid bad=g; bad.spacing.x()=0; dualContour(bad,sphere); } catch(const std::invalid_argument&) { rejected=true; }
    check(rejected,"Invalid grid accepted");
    rejected=false;
    try { dualContour(g,sphere,[](const Vec3&)->Vec3{return Vec3::Zero();}); } catch(const std::invalid_argument&) { rejected=true; }
    check(rejected,"Zero normals accepted");
    auto output=std::filesystem::path("test-sphere.obj");
    writeObj(mesh,output); std::ifstream in(output); std::string first; in>>first;
    check(first=="v","Invalid OBJ output"); in.close(); std::filesystem::remove(output);
    std::cout<<"Sphere topology/winding, non-cubic grid, plane QEF, empty field, input validation and OBJ passed\n";
 } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}

