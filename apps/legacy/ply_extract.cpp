#include "PlyDistanceField.h"
#include "PostTopOpt.h"
#include "ObjWriter.h"
#include <chrono>
#include <iostream>
#include <string>

int main(int argc,char** argv) {
    try {
        std::filesystem::path input,output;
        std::string methodName="dc";
        int resolution=64;
        for(int i=1;i<argc;++i) {
            std::string arg=argv[i];
            if(arg=="--help") {
                std::cout<<"Usage: PostTopOptPlyExtract input.ply [--output file.obj] [--cells 64] [--method dc|emilk|mc]\n"
                    "Closed inputs only; --cells counts cells across the longest unpadded side.\n";
                return 0;
            }
            if(arg=="--output" || arg=="--cells" || arg=="--method") {
                if(++i==argc) throw std::invalid_argument("Missing option value");
                if(arg=="--output") output=argv[i];
                else if(arg=="--method") methodName=argv[i];
                else { std::size_t used=0; resolution=std::stoi(argv[i],&used);
                    if(used!=std::string(argv[i]).size()) throw std::invalid_argument("Invalid cell count"); }
            } else if(arg.starts_with("--") || !input.empty()) throw std::invalid_argument("Unexpected argument: "+arg);
            else input=arg;
        }
        if(input.empty()) throw std::invalid_argument("Input PLY required; use --help");
        if(methodName!="dc" && methodName!="emilk" && methodName!="mc")
            throw std::invalid_argument("This sequential VTK example supports dc, emilk and mc");
        const auto method=geometry::parseMethod(methodName);
        if(output.empty()) output=std::filesystem::path(POSTTOPOPT_RESULTS_DIR)/"stanford"/
            (input.stem().string()+"_"+methodName+"_"+std::to_string(resolution)+".obj");
        if(std::filesystem::exists(output)) throw std::invalid_argument("Output exists; choose a new --output path");
        const auto begin=std::chrono::steady_clock::now();
        auto model=geometry::makePlyDistanceField(input,resolution);
        geometry::ContourOptions options;
        options.gradientStep=model.grid.spacing.minCoeff()*1e-3;
        std::cout<<"source_triangles="<<model.source.sourcePolygons
            <<"\norigin="<<model.grid.origin.transpose()<<"\nspacing="<<model.grid.spacing.transpose()
            <<"\ncells="<<model.grid.cells[0]<<','<<model.grid.cells[1]<<','<<model.grid.cells[2]
            <<"\nregularization="<<options.regularization<<", constrained="<<options.constrainToCell
            <<", bisection="<<options.bisection<<", gradientStep="<<options.gradientStep<<std::endl;
        const auto extractionBegin=std::chrono::steady_clock::now();
        auto mesh=geometry::extractSurface(method,model.grid,model.field,{},options);
        const auto end=std::chrono::steady_clock::now();
        if(mesh.vertices.empty() || (mesh.quads.empty() && mesh.triangles.empty()))
            throw std::runtime_error("Extraction produced an empty surface");
        if(!output.parent_path().empty()) std::filesystem::create_directories(output.parent_path());
        geometry::writeObj(mesh,output);
        std::cout<<"method="<<methodName<<", vertices="<<mesh.vertices.size()<<", quads="<<mesh.quads.size()
            <<", triangles="<<mesh.triangles.size()<<", extraction_ms="
            <<std::chrono::duration<double,std::milli>(end-extractionBegin).count()<<", preparation_ms="
            <<std::chrono::duration<double,std::milli>(extractionBegin-begin).count()<<" -> "<<output.string()<<'\n';
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
