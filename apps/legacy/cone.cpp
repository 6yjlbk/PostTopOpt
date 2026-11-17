#include "ExampleRunner.h"
#include <algorithm>
#include <cmath>

int main(int argc, char** argv) {
    // Same cone as gradient_comparison.cpp: R=1, H=2, base z=0.
    auto grid=examples::grid(1.5,80);
    grid.origin.z()=-0.5;
    const auto field=[](const geometry::Vec3& p) {
        return std::max(std::hypot(p.x(),p.y())+0.5*p.z()-1.0,-p.z());
    };
#ifdef POSTTOPOPT_CONE_ANALYTIC
    const auto gradient=[](const geometry::Vec3& p)->geometry::Vec3 {
        const double r=std::hypot(p.x(),p.y());
        const double side=r+0.5*p.z()-1.0;
        // Nondifferentiable tie: select the base branch.
        if (-p.z()>=side) return geometry::Vec3(0,0,-1);
        // Axis: symmetric subgradient, not a unique classical derivative.
        if (r==0) return geometry::Vec3(0,0,0.5);
        return geometry::Vec3(p.x()/r,p.y()/r,0.5);
    };
    return examples::run(argc,argv,"cone_analytic.obj",grid,field,gradient);
#else
    return examples::run(argc,argv,"cone_numeric.obj",grid,field);
#endif
}
