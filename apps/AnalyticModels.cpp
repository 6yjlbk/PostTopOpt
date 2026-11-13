#include "AnalyticModels.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace application {
using namespace geometry;
CellScalarVolume generateAnalytic(const std::string& name,int n,int complexity){
    if(n<16 || n>500 || complexity<1 || complexity>12) throw std::invalid_argument("cells:16..500, complexity:1..12");
    Grid g;g.origin=Vec3::Constant(-1.);g.spacing=Vec3::Constant(2./n);g.cells={n,n,n};
    // Smooth synthetic fractions, not exact volume integrals of a binary solid.
    return sampleCells(g,[&](const Vec3& p){
        const double k=complexity*3.14159265358979323846;
        const double envelope=.85-p.cwiseAbs().maxCoeff();
        double margin=0;
        if(name=="gyroid"){
            const double t=std::sin(k*p.x())*std::cos(k*p.y())+std::sin(k*p.y())*std::cos(k*p.z())+std::sin(k*p.z())*std::cos(k*p.x());
            margin=(.35-std::abs(t))/k;
        }else if(name=="channels"){
            const double x=std::sin(k*p.x()),y=std::sin(k*p.y()),z=std::sin(k*p.z());
            margin=(std::min({std::hypot(x,y),std::hypot(y,z),std::hypot(z,x)})-.38)/k;
        }else if(name=="lattice"){
            const double x=std::sin(k*p.x()),y=std::sin(k*p.y()),z=std::sin(k*p.z());
            margin=(.30-std::min({std::hypot(x,y),std::hypot(y,z),std::hypot(z,x)}))/k;
        }else if(name=="shells") margin=(.22-std::abs(std::sin(k*p.norm())))/k;
        else if(name=="sphere") margin=.65-p.norm();
        else if(name=="torus") margin=.18-std::hypot(std::hypot(p.x(),p.y())-.48,p.z());
        else if(name=="sphere-union") margin=std::max(.45-(p-Vec3(-.25,0,0)).norm(),.45-(p-Vec3(.25,0,0)).norm());
        else if(name=="drilled-box") margin=std::min(.65-p.cwiseAbs().maxCoeff(),std::min({std::hypot(p.x(),p.y()),std::hypot(p.y(),p.z()),std::hypot(p.z(),p.x())})-.22);
        else if(name=="cone") margin=std::min({.5-std::hypot(p.x(),p.y())-.5*p.z(),p.z()+.5,.5-p.z()});
        else throw std::invalid_argument("Unknown analytic model: "+name);
        return std::clamp(.5+std::min(margin,envelope)/.04,0.,1.);
    });
}
}
