#include "DualContouring.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
using namespace geometry;
void compare(const Mesh& a,const Mesh& b){
    if(a.quads!=b.quads || a.vertices.size()!=b.vertices.size())throw std::runtime_error("Optimization changed connectivity");
    for(std::size_t i=0;i<a.vertices.size();++i)
        if((a.vertices[i]-b.vertices[i]).norm()>1e-8*std::max(1.,b.vertices[i].norm()))throw std::runtime_error("Optimization changed coordinates");
}
int main(){try{
    Grid g;g.origin=Vec3(-.7,-.8,-.9);g.spacing=Vec3(.17,.21,.23);g.cells={8,7,6};
    std::mt19937 engine(4261);std::uniform_real_distribution<double> distribution(-1.,1.);
    std::size_t checked=0,vertices=0;
    for(int trial=0;trial<12;++trial){
        Vec3 linear(distribution(engine),distribution(engine),distribution(engine));
        const double c=distribution(engine)*.2;
        Field f=[=](const Vec3& p){return linear.dot(p)+.3*p.squaredNorm()+.1*std::sin(3*p.x())-c;};
        Gradient gradient=[=](const Vec3& p)->Vec3 {return linear+.6*p+Vec3(.3*std::cos(3*p.x()),0,0);};
        for(bool constrained:{false,true})for(bool regularized:{false,true})for(bool bisect:{false,true}){
            ContourOptions original;original.constrainToCell=constrained;original.regularizationEnabled=regularized;original.bisection=bisect;
            original.optimizeQef=false;original.cacheHermite=false;
            auto optimized=original;optimized.optimizeQef=true;optimized.cacheHermite=true;
            auto a=dualContour(g,f,gradient,original),b=dualContour(g,f,gradient,optimized);compare(a,b);vertices+=b.vertices.size();++checked;
            compare(dualContourEmilk(g,f,gradient,original),dualContourEmilk(g,f,gradient,optimized));++checked;
        }
    }
    // Rank-deficient plane, exact grid-node roots, endpoint interpolation.
    for(double iso:{0.,.15}){
        ContourOptions a;a.isovalue=iso;a.regularization=0;a.optimizeQef=false;a.cacheHermite=false;
        auto b=a;b.optimizeQef=true;b.cacheHermite=true;
        compare(dualContour(g,[](const Vec3& p){return p.x();},{},a),dualContour(g,[](const Vec3& p){return p.x();},{},b));
    }
    if(!vertices)throw std::runtime_error("No active cells checked");
    std::cout<<checked<<" optimized/reference comparisons passed, "<<vertices<<" DC vertices checked\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
