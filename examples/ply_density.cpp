#include "PlyDistanceField.h"
#include "PostTopOpt.h"
#include "ObjWriter.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
namespace {
using Clock=std::chrono::steady_clock;
double milliseconds(Clock::time_point start){return std::chrono::duration<double,std::milli>(Clock::now()-start).count();}
geometry::Backend backend(const std::string& s){
    if(s=="dc")return geometry::Backend::Dc;
    if(s=="emilk")return geometry::Backend::Emilk;
    if(s=="mc")return geometry::Backend::Mc;
    if(s=="cgal-dc")return geometry::Backend::CgalDc;
    if(s=="paraview-amr")return geometry::Backend::ParaViewAmr;
    throw std::invalid_argument("Unknown method: "+s);
}
}
int main(int argc,char** argv){try{
    std::filesystem::path input, directory;
    std::string method="dc";
    int resolution=32;
    double widthCells=2.;
    geometry::ExtractionSettings settings;
    for(int i=1;i<argc;++i){
        std::string arg=argv[i];
        if(arg=="--help"){
            std::cout<<"Usage: PostTopOptPlyDensity input.ply --output-dir directory [--cells 32] [--width-cells 2] [--method all|dc|emilk|mc|cgal-dc|paraview-amr] [--lambda 0.001] [--regularization 0|1] [--constrained 0|1] [--bisection 0|1]\n";return 0;
        }
        if(!arg.starts_with("--")){if(!input.empty())throw std::invalid_argument("Multiple input files");input=arg;continue;}
        if(++i==argc)throw std::invalid_argument("Missing option value");
        std::string value=argv[i];std::size_t used=0;
        if(arg=="--output-dir")directory=value;
        else if(arg=="--method")method=value;
        else if(arg=="--cells"){resolution=std::stoi(value,&used);if(used!=value.size())throw std::invalid_argument("Invalid cells");}
        else if(arg=="--width-cells"){widthCells=std::stod(value,&used);if(used!=value.size())throw std::invalid_argument("Invalid width");}
        else if(arg=="--lambda"){settings.dc.regularization=std::stod(value,&used);if(used!=value.size() || !std::isfinite(settings.dc.regularization) || settings.dc.regularization<0)throw std::invalid_argument("Invalid lambda");}
        else {
            if(value!="0" && value!="1")throw std::invalid_argument("Boolean requires 0 or 1");
            bool flag=value=="1";
            if(arg=="--regularization")settings.dc.regularizationEnabled=flag;
            else if(arg=="--constrained"){settings.dc.constrainToCell=flag;settings.cgal.constrainToCell=flag;}
            else if(arg=="--bisection"){settings.dc.bisection=flag;settings.cgal.bisection=flag;}
            else throw std::invalid_argument("Unknown option: "+arg);
        }
    }
    if(input.empty() || directory.empty())throw std::invalid_argument("Input and --output-dir required");
    // The SDF adapter reserves 3 cells; cap the transition width to retain an empty boundary.
    if(!std::isfinite(widthCells) || widthCells<=0 || widthCells>4)throw std::invalid_argument("width-cells must be in (0,4]");
    std::vector<std::string> methods=method=="all"?std::vector<std::string>{"dc","emilk","mc","cgal-dc","paraview-amr"}:std::vector<std::string>{method};
    for(const auto& name:methods){
        backend(name);
        if(name=="cgal-dc" && !geometry::cgalAvailable())throw std::invalid_argument("CGAL not enabled");
        if(name=="paraview-amr" && !geometry::paraViewAvailable())throw std::invalid_argument("ParaView not enabled");
        if(std::filesystem::exists(directory/(name+".obj")))throw std::invalid_argument("Output exists; choose a new directory");
    }
    if(std::filesystem::exists(directory/"report.csv"))throw std::invalid_argument("Report exists; choose a new directory");
    const auto start=Clock::now();
    auto model=geometry::makePlyDistanceField(input,resolution);
    const double preparationMs=milliseconds(start);
    const double h=model.grid.spacing.minCoeff(), width=h*widthCells;
    const auto samplingStart=Clock::now();
    auto volume=geometry::sampleCells(model.grid,[&](const geometry::Vec3& p){return std::clamp(.5-model.field(p)/width,0.,1.);});
    const double samplingMs=milliseconds(samplingStart);
    // Release the expensive source geometry before extraction; all methods use this same array.
    model.field={};
    settings.dc.gradientStep=h*1e-3; settings.cgal.gradientStep=h*1e-3;
    std::filesystem::create_directories(directory);
    std::ofstream report(directory/"report.csv");
    if(!report)throw std::runtime_error("Cannot create report");
    report<<std::setprecision(17);
    report<<"method,source,cells_requested,nx,ny,nz,ox,oy,oz,h,width_cells,width,threshold,lambda,regularization,constrained,bisection,gradient_step,preparation_ms,sampling_ms,extraction_ms,export_ms,vertices,quads,triangles\n";
    for(const auto& name:methods){
        const auto extractionStart=Clock::now();
        auto mesh=geometry::extract(backend(name),volume,.5,settings);
        const double extractionMs=milliseconds(extractionStart);
        if(mesh.vertices.empty() || (mesh.quads.empty() && mesh.triangles.empty()))throw std::runtime_error("Empty extraction: "+name);
        const auto exportStart=Clock::now();geometry::writeObj(mesh,directory/(name+".obj"));
        const double exportMs=milliseconds(exportStart);
        std::cout<<name<<": "<<mesh.vertices.size()<<" vertices, "<<mesh.quads.size()<<" quads, "<<mesh.triangles.size()<<" triangles, extraction_ms="<<extractionMs<<'\n';
        report<<name<<','<<std::quoted(input.string())<<','<<resolution<<','<<model.grid.cells[0]<<','<<model.grid.cells[1]<<','<<model.grid.cells[2]<<','
            <<model.grid.origin.x()<<','<<model.grid.origin.y()<<','<<model.grid.origin.z()<<','<<h<<','<<widthCells<<','<<width<<",0.5,"<<settings.dc.regularization<<','
            <<settings.dc.regularizationEnabled<<','<<settings.dc.constrainToCell<<','<<settings.dc.bisection<<','<<settings.dc.gradientStep<<','
            <<preparationMs<<','<<samplingMs<<','<<extractionMs<<','<<exportMs<<','<<mesh.vertices.size()<<','<<mesh.quads.size()<<','<<mesh.triangles.size()<<'\n';
        report.flush();if(!report)throw std::runtime_error("Cannot write report");
    }
    std::cout<<"preparation_ms="<<preparationMs<<", sampling_ms="<<samplingMs<<", width="<<width<<", report="<<(directory/"report.csv").string()<<'\n';
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
