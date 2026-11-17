#include "PostTopOpt.h"
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}
int main(){try{
    geometry::Grid grid;grid.origin=geometry::Vec3(-1,-2,-3);grid.spacing=geometry::Vec3(.2,.3,.4);grid.cells={4,3,2};
    auto field=[](const geometry::Vec3& p){return 2*p.x()+3*p.y()-p.z();};
    auto volume=std::make_shared<const geometry::CellScalarVolume>(geometry::sampleCells(grid,field));
    auto reconstructed=geometry::linearCellField(volume);
    geometry::Vec3 p=grid.origin+grid.spacing.cwiseProduct(geometry::Vec3(1.2,1.3,.8));
    require(std::abs(reconstructed(p)-field(p))<1e-12,"Linear reproduction failed");
    require(geometry::cellIndex(grid,1,0,0)==6,"Wrong memory layout");
    auto rho=geometry::densityField(volume,.5);
    require(std::abs(rho(p)-(.5-field(p)))<1e-12,"Wrong density sign");
    bool rejected=false;try{auto bad=*volume;bad.values.pop_back();bad.validate();}catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"Bad array size accepted");
    geometry::ExtractionSettings options;
    require(options.dc.regularization==.001 && options.dc.regularizationEnabled && options.dc.constrainToCell && options.dc.bisection,"DC baseline changed");
    geometry::Grid sphere; sphere.origin=geometry::Vec3::Constant(-1.5);sphere.spacing=geometry::Vec3::Constant(.25);sphere.cells={12,12,12};
    auto mesh=geometry::extract(geometry::Backend::Dc,sphere,[](const geometry::Vec3& q){return q.norm()-1.;});
    require(!mesh.quads.empty(),"Analytic API regression");
    options.dc.regularizationEnabled=false;
    auto plane=geometry::extract(geometry::Backend::Dc,sphere,[](const geometry::Vec3& q){return q.x()-.1;},{},options);
    for(const auto& q:plane.vertices) require(std::abs(q.x()-.1)<1e-3,"Regularization toggle failed");
    if(!geometry::paraViewAvailable()){
        rejected=false;try{geometry::paraViewDualContour(*volume);}catch(const std::runtime_error&){rejected=true;}
        require(rejected,"Disabled backend did not report error");
    }
    if(geometry::paraViewAvailable()) {
        auto fractions=geometry::sampleCells(sphere,[](const geometry::Vec3& q){return q.norm()<1. ? 1. : 0.;});
        auto result=geometry::extract(geometry::Backend::ParaViewAmr,fractions,.5);
        require(!result.vertices.empty() && !result.triangles.empty(),"ParaView extraction failed");
        for(const auto& p:result.vertices) require(p.allFinite(),"Non-finite ParaView point");
        for(const auto& t:result.triangles) for(auto id:t) require(id<result.vertices.size(),"Invalid ParaView face");
    }
    std::cout<<"Architecture API checks passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}

