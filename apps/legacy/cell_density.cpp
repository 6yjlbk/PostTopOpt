#include <algorithm>
#include "PostTopOpt.h"
#include "ObjWriter.h"
#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
int main(int argc,char** argv) {
    try {
        geometry::Backend backend=geometry::Backend::Dc;
        geometry::ExtractionSettings settings;
        std::filesystem::path output="cell-density.obj";
        for(int i=1;i<argc;++i) {
            std::string arg=argv[i];
            if(arg=="--help") {std::cout<<"Usage: PostTopOptCellDensity [--method dc|emilk|mc|cgal-dc|paraview-amr] [--output path] [--lambda value] [--regularization 0|1] [--constrained 0|1] [--bisection 0|1]\n";return 0;}
            if(++i==argc) throw std::invalid_argument("Missing option value");
            std::string value=argv[i];
            if(arg=="--method") {
                if(value=="dc") backend=geometry::Backend::Dc;
                else if(value=="emilk") backend=geometry::Backend::Emilk;
                else if(value=="mc") backend=geometry::Backend::Mc;
                else if(value=="cgal-dc") backend=geometry::Backend::CgalDc;
                else if(value=="paraview-amr") backend=geometry::Backend::ParaViewAmr;
                else throw std::invalid_argument("Unknown backend");
            } else if(arg=="--output") output=value;
            else if(arg=="--lambda") {std::size_t used; settings.dc.regularization=std::stod(value,&used); if(used!=value.size() || !std::isfinite(settings.dc.regularization) || settings.dc.regularization<0) throw std::invalid_argument("Invalid lambda");}
            else {
                if(value!="0" && value!="1") throw std::invalid_argument("Boolean options accept 0 or 1");
                bool enabled=value=="1";
                if(arg=="--regularization") settings.dc.regularizationEnabled=enabled;
                else if(arg=="--constrained") {settings.dc.constrainToCell=enabled;settings.cgal.constrainToCell=enabled;}
                else if(arg=="--bisection") {settings.dc.bisection=enabled;settings.cgal.bisection=enabled;}
                else throw std::invalid_argument("Unknown option");
            }
        }
        geometry::Grid grid; grid.origin=geometry::Vec3::Constant(-1.5);
        grid.spacing=geometry::Vec3::Constant(0.1);grid.cells={30,30,30};
        // Synthetic cell fractions, not exact integrated material volumes.
        auto volume=geometry::sampleCells(grid,[](const geometry::Vec3& p){return std::clamp(0.5+(1.-p.norm())/0.2,0.,1.);});
        const auto mesh=geometry::extract(backend,volume,0.5,settings);
        if(!output.parent_path().empty()) std::filesystem::create_directories(output.parent_path());
        if(std::filesystem::exists(output)) throw std::invalid_argument("Output already exists");
        geometry::writeObj(mesh,output);
        std::cout<<mesh.vertices.size()<<" vertices, "<<mesh.quads.size()<<" quads, "<<mesh.triangles.size()<<" triangles -> "<<output.string()<<'\n';
    } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
