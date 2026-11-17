#include "VolumeIO.h"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
namespace application {
using namespace geometry;
CellScalarVolume loadPtov(const std::filesystem::path& p){
    std::ifstream f(p);std::string magic;f>>magic;
    if(magic!="PTOV1")throw std::runtime_error("Invalid volume format");
    CellScalarVolume v;
    for(int d=0;d<3;++d)f>>v.grid.origin[d];
    for(int d=0;d<3;++d)f>>v.grid.spacing[d];
    for(int& n:v.grid.cells)f>>n;
    if(!f)throw std::runtime_error("Invalid volume header");validateGrid(v.grid);
    const auto count=cellCount(v.grid);if(count>500ULL*256*256)throw std::runtime_error("Benchmark volume limit exceeded");
    v.values.resize(count);for(double& value:v.values)if(!(f>>value))throw std::runtime_error("Truncated volume");
    f>>std::ws;if(!f.eof())throw std::runtime_error("Unexpected trailing volume data");
    v.validate();for(double value:v.values)if(value<0 || value>1)throw std::runtime_error("Fractions must be in [0,1]");return v;
}

// CSV schema: version, grid geometry, then explicitly indexed cell values.
void saveVolume(const CellScalarVolume& v,const std::filesystem::path& path) {
    v.validate();
    if(std::filesystem::exists(path)) throw std::runtime_error("Volume already exists");
    if(!path.parent_path().empty())std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path);out<<std::setprecision(17);
    out<<"PostTopOptVolume,1\n";
    out<<"cells,"<<v.grid.cells[0]<<','<<v.grid.cells[1]<<','<<v.grid.cells[2]<<'\n';
    out<<"origin,"<<v.grid.origin.x()<<','<<v.grid.origin.y()<<','<<v.grid.origin.z()<<'\n';
    out<<"spacing,"<<v.grid.spacing.x()<<','<<v.grid.spacing.y()<<','<<v.grid.spacing.z()<<"\nx,y,z,rho\n";
    for(int x=0;x<v.grid.cells[0];++x)for(int y=0;y<v.grid.cells[1];++y)for(int z=0;z<v.grid.cells[2];++z)
        out<<x<<','<<y<<','<<z<<','<<v.values[cellIndex(v.grid,x,y,z)]<<'\n';
    if(!out)throw std::runtime_error("Cannot write density CSV");
}
CellScalarVolume loadVolume(const std::filesystem::path& path) {
    if(path.extension()==".ptov")return loadPtov(path);
    std::ifstream in(path);std::string line;
    auto row=[&](){if(!std::getline(in,line))throw std::runtime_error("Truncated density CSV");if(!line.empty() && line.back()=='\r')line.pop_back();std::replace(line.begin(),line.end(),',',' ');return std::istringstream(line);};
    {auto r=row();std::string magic;int version=0;r>>magic>>version;if(magic!="PostTopOptVolume" || version!=1)throw std::runtime_error("Invalid density CSV version");}
    CellScalarVolume v;
    auto header=[&](const char* expected,auto& values){auto r=row();std::string key;r>>key;for(int d=0;d<3;++d)r>>values[d];if(!r)throw std::runtime_error("Invalid density CSV number");r>>std::ws;if(key!=expected || !r.eof())throw std::runtime_error("Invalid density CSV header");};
    header("cells",v.grid.cells);header("origin",v.grid.origin);header("spacing",v.grid.spacing);
    validateGrid(v.grid);if(cellCount(v.grid)>500ULL*256*256)throw std::runtime_error("Volume limit: 500^3 cells");
    {auto r=row();std::string a,b,c,d;r>>a>>b>>c>>d;if(a!="x" || b!="y" || c!="z" || d!="rho")throw std::runtime_error("Invalid density CSV columns");}
    v.values.resize(cellCount(v.grid));
    for(int x=0;x<v.grid.cells[0];++x)for(int y=0;y<v.grid.cells[1];++y)for(int z=0;z<v.grid.cells[2];++z){auto r=row();int i=-1,j=-1,k=-1;double rho=-1;r>>i>>j>>k>>rho;if(!r)throw std::runtime_error("Invalid density CSV number");r>>std::ws;if(!r.eof() || i!=x || j!=y || k!=z || rho<0 || rho>1)throw std::runtime_error("Invalid density CSV cell");v.values[cellIndex(v.grid,x,y,z)]=rho;}
    in>>std::ws;if(!in.eof())throw std::runtime_error("Trailing density CSV data");v.validate();return v;
}
}
