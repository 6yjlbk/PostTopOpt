#include "CgalExtraction.h"
#include <Eigen/Geometry>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace geometry;
void check(bool v,const char* m) {if(!v) throw std::runtime_error(m);}
int main() {
 try {
    Grid g;g.origin=Vec3::Zero();g.spacing=Vec3(.2,.3,.4);g.cells={4,3,2};
    auto field=[](const Vec3& p){return p.x();};
    if(!cgalAvailable()) {
        try {cgalDualContour(g,field);} catch(const std::runtime_error& e) {
            check(std::string(e.what()).find("POSTTOPOPT_WITH_CGAL=ON")!=std::string::npos,"Missing disabled message");
            std::cout<<"Disabled CGAL reports the required build option\n";return 0;
        }
        throw std::runtime_error("Disabled backend silently succeeded");
    }
    CgalContourOptions o;o.isovalue=.35;
    auto m=cgalDualContour(g,field,{},o);
    check(m.vertices.size()==6 && m.quads.size()==2,"Grid cell/node conversion mismatch");
    for(const auto& p:m.vertices) check(p.allFinite() && std::abs(p.x()-.35)<1e-8,"Wrong plane");
    for(const auto& q:m.quads) {
        for(auto i:q) check(i<m.vertices.size(),"Bad index");
        check((m.vertices[q[1]]-m.vertices[q[0]]).cross(m.vertices[q[2]]-m.vertices[q[0]]).x()>0,"Wrong winding");
    }
    o.constrainToCell=false;
    auto analytical=cgalDualContour(g,field,[](const Vec3&)->Vec3{return Vec3(1,0,0);},o);
    for(const auto& p:analytical.vertices) check(std::abs(p.x()-.35)<1e-8,"Free/analytic plane");
    check(cgalDualContour(g,[](const Vec3&){return 1.;}).vertices.empty(),"Constant field not empty");
    bool rejected=false;
    try {auto bad=g;bad.spacing.x()=0;cgalDualContour(bad,field);}catch(const std::invalid_argument&){rejected=true;}
    check(rejected,"Invalid grid accepted");
    rejected=false;
    try {cgalDualContour(g,field,[](const Vec3&)->Vec3{return Vec3::Zero();},o);}catch(const std::invalid_argument&){rejected=true;}
    check(rejected,"Invalid normal accepted");
    std::cout<<"CGAL plane, anisotropic grid, winding, gradients, free/constrained and validation passed\n";
 } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
