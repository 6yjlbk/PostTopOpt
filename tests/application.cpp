#include "AnalyticModels.h"
#include "VolumeIO.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
int main(int argc,char** argv) {try {
    if(argc!=2)throw std::runtime_error("Test directory required");
    const std::filesystem::path directory=argv[1];std::filesystem::create_directories(directory);
    for(const char* name:{"sphere","torus","sphere-union","drilled-box","cone","gyroid","channels","lattice","shells"}) {
        auto volume=application::generateAnalytic(name,16,2);volume.validate();
        for(double rho:volume.values)if(!std::isfinite(rho) || rho<0 || rho>1)throw std::runtime_error("Invalid synthetic density");
    }
    auto volume=application::generateAnalytic("sphere",16,2);
    const auto path=directory/"roundtrip.csv";
    // This path belongs exclusively to this test.
    std::filesystem::remove(path);
    application::saveVolume(volume,path);const auto loaded=application::loadVolume(path);
    if(volume.values!=loaded.values || volume.grid.cells!=loaded.grid.cells ||
       !volume.grid.origin.isApprox(loaded.grid.origin,0) || !volume.grid.spacing.isApprox(loaded.grid.spacing,0))throw std::runtime_error("CSV roundtrip changed volume");
    {std::ofstream out(directory/"invalid.csv");out<<"PostTopOptVolume,1\ncells,2,2\n";}
    bool rejected=false;try{application::loadVolume(directory/"invalid.csv");}catch(const std::exception&){rejected=true;}
    if(!rejected)throw std::runtime_error("Malformed CSV accepted");
    std::cout<<"Application data checks passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
