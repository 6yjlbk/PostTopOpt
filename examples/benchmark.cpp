#include "PostTopOpt.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>
namespace {
using namespace geometry;
using Clock=std::chrono::steady_clock;
double ms(Clock::time_point t){return std::chrono::duration<double,std::milli>(Clock::now()-t).count();}
void save(const CellScalarVolume& v,const std::filesystem::path& p){
    std::ofstream f(p); if(!f)throw std::runtime_error("Cannot save volume");
    f<<"PTOV1\n"<<std::setprecision(17)<<v.grid.origin.transpose()<<'\n'<<v.grid.spacing.transpose()<<'\n';
    for(int n:v.grid.cells)f<<n<<' ';f<<'\n';
    for(double value:v.values)f<<value<<'\n';
    if(!f)throw std::runtime_error("Volume write failed");
}
CellScalarVolume load(const std::filesystem::path& p){
    std::ifstream f(p);std::string magic;f>>magic;
    if(magic!="PTOV1")throw std::runtime_error("Invalid volume format");
    CellScalarVolume v;
    for(int d=0;d<3;++d)f>>v.grid.origin[d];
    for(int d=0;d<3;++d)f>>v.grid.spacing[d];
    for(int& n:v.grid.cells)f>>n;
    if(!f)throw std::runtime_error("Invalid volume header");validateGrid(v.grid);
    const auto count=cellCount(v.grid);if(count>256ULL*256*256)throw std::runtime_error("Benchmark volume limit exceeded");
    v.values.resize(count);for(double& value:v.values)if(!(f>>value))throw std::runtime_error("Truncated volume");
    f>>std::ws;if(!f.eof())throw std::runtime_error("Unexpected trailing volume data");
    v.validate();for(double value:v.values)if(value<0 || value>1)throw std::runtime_error("Fractions must be in [0,1]");return v;
}
CellScalarVolume generate(const std::string& name,int n,int complexity){
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
        else throw std::invalid_argument("Unknown case");
        return std::clamp(.5+std::min(margin,envelope)/.04,0.,1.);
    });
}
struct GeometrySummary{std::size_t vertices,quads,triangles;};
GeometrySummary summary(const Mesh& mesh){
    if(mesh.vertices.empty() || (mesh.quads.empty() && mesh.triangles.empty()))throw std::runtime_error("Empty benchmark output");
    for(const auto& p:mesh.vertices)if(!p.allFinite())throw std::runtime_error("Non-finite output");
    return {mesh.vertices.size(),mesh.quads.size(),mesh.triangles.size()};
}
double median(std::vector<double> v){std::sort(v.begin(),v.end());return (v[(v.size()-1)/2]+v[v.size()/2])/2.;}
}
int main(int argc,char** argv){try{
    std::string task="run";std::filesystem::path directory;int n=48,complexity=3,repeats=5,cache=1,qef=1;
    for(int i=1;i<argc;++i){std::string arg=argv[i];
        if(arg=="--help"){std::cout<<"PostTopOptBenchmark --task generate|run --bank directory [--cells 48] [--complexity 3] [--repeats 5] [--cache 0|1] [--qef 0|1]\n";return 0;}
        if(++i==argc)throw std::invalid_argument("Missing argument");std::string value=argv[i];std::size_t used=0;
        if(arg=="--task")task=value;else if(arg=="--bank")directory=value;
        else{int parsed=std::stoi(value,&used);if(used!=value.size())throw std::invalid_argument("Invalid integer");
            if(arg=="--cells")n=parsed;else if(arg=="--complexity")complexity=parsed;else if(arg=="--repeats")repeats=parsed;else if(arg=="--cache")cache=parsed;else if(arg=="--qef")qef=parsed;else throw std::invalid_argument("Unknown argument");}
    }
    if((qef!=0 && qef!=1) || (cache!=0 && cache!=1))throw std::invalid_argument("cache and qef:0|1");
    if(directory.empty())throw std::invalid_argument("--bank required");
    if(task=="generate"){
        if(n<16 || n>256 || complexity<1 || complexity>12)throw std::invalid_argument("cells:16..256 complexity:1..12");
        const std::vector<std::string> cases={"gyroid","channels","lattice","shells"};
        for(const auto& name:cases)if(std::filesystem::exists(directory/(name+".ptov")))throw std::invalid_argument("Bank already exists");
        std::filesystem::create_directories(directory);
        std::ofstream manifest(directory/"manifest.txt");manifest<<"generator=1\ncells="<<n<<"\ncomplexity="<<complexity<<"\nramp_width=0.04\norigin=-1,-1,-1\nindex=z-fastest\n";
        for(const auto& name:cases){save(generate(name,n,complexity),directory/(name+".ptov"));std::cout<<"generated "<<name<<'\n';}
        return 0;
    }
    if(task!="run" || repeats<2 || repeats>50)throw std::invalid_argument("task run requires repeats:2..50");
    if(!geometry::paraViewAvailable())throw std::runtime_error("ParaView backend required");
    std::vector<std::filesystem::path> inputs;
    for(const auto& entry:std::filesystem::directory_iterator(directory))if(entry.is_regular_file() && entry.path().extension()==".ptov")inputs.push_back(entry.path());
    std::sort(inputs.begin(),inputs.end());if(inputs.empty())throw std::runtime_error("Empty bank");
    if(std::filesystem::exists(directory/"runs.csv") || std::filesystem::exists(directory/"summary.csv"))throw std::runtime_error("Reports exist; move them before repeating");
    std::ofstream config(directory/"run-settings.txt");config<<"cache="<<cache<<"\nqef="<<qef<<"\nrepeats="<<repeats<<"\nthreshold=0.5\nlambda=0.001\nregularization=1\nconstrained=1\nbisection=1\n";
    std::ofstream raw(directory/"runs.csv"),out(directory/"summary.csv");if(!raw || !out)throw std::runtime_error("Cannot create reports");
    raw<<"case,method,repeat,total_ms,vertices,quads,triangles\n";
    out<<"case,method,cells,total_median_ms,total_min_ms,total_max_ms,profile_total_ms,prepare_ms,node_sampling_ms,hermite_ms,qef_ms,connectivity_ms,filter_ms,conversion_ms,active_cells,edge_visits,unique_edges,field_queries,qef_interior,qef_boundary,hermite_evaluations,hermite_cache_hits,vertices,quads,triangles\n";
    raw<<std::setprecision(10);out<<std::setprecision(10);
    for(const auto& path:inputs){
        auto volume=load(path); // Disk I/O is outside timing.
        const auto caseName=path.stem().string();
        GeometrySummary counts{};std::vector<double> times[2];
        auto run=[&](int method,const ExtractionSettings& settings){return extract(method==0?Backend::Dc:Backend::ParaViewAmr,volume,.5,settings);};
        ExtractionSettings settings;settings.dc.cacheHermite=cache!=0;settings.dc.optimizeQef=qef!=0;settings.dc.gradientStep=volume.grid.spacing.minCoeff()*1e-3;
        // A/B geometry check with coordinate tolerance outside timing. This also covers constrained QEF.
        {
            auto cached=settings,reference=settings;
            cached.dc.cacheHermite=true; reference.dc.cacheHermite=false;reference.dc.optimizeQef=false;
            DcProfile pc,pr; cached.dc.profile=&pc; reference.dc.profile=&pr;
            auto a=run(0,cached),b=run(0,reference);
            if(a.quads!=b.quads || a.triangles!=b.triangles || a.vertices.size()!=b.vertices.size())
                throw std::runtime_error("Hermite cache changed mesh connectivity");
            for(std::size_t i=0;i<a.vertices.size();++i)
                if((a.vertices[i]-b.vertices[i]).norm()>1e-10*std::max(1.,volume.grid.spacing.maxCoeff()))throw std::runtime_error("Hermite cache changed vertex coordinates");
            if(pc.hermiteEvaluations+pc.hermiteCacheHits!=pc.edgeVisits || pc.fieldQueries>pr.fieldQueries)
                throw std::runtime_error("Inconsistent Hermite cache counters");
        }
        // Warm up both backends, then alternate order to reduce order bias.
        for(int method=0;method<2;++method)summary(run(method,settings));
        GeometrySummary expected[2]{};
        for(int rep=0;rep<repeats;++rep)for(int order=0;order<2;++order){int method=(rep+order)%2;
            const auto begin=Clock::now();auto mesh=run(method,settings);const double total=ms(begin);
            counts=summary(mesh); // Validation is outside timing; allocation/conversion are included.
            if(rep && (counts.vertices!=expected[method].vertices || counts.quads!=expected[method].quads || counts.triangles!=expected[method].triangles))throw std::runtime_error("Non-deterministic output counts");
            expected[method]=counts;times[method].push_back(total);
            raw<<caseName<<','<<(method==0?"dc":"paraview-amr")<<','<<rep<<','<<total<<','<<counts.vertices<<','<<counts.quads<<','<<counts.triangles<<'\n';
        }
        for(int method=0;method<2;++method){
            DcProfile dc;ParaViewProfile pv;auto profiled=settings;
            if(method==0)profiled.dc.profile=&dc;else profiled.paraview.profile=&pv;
            const auto begin=Clock::now();auto mesh=run(method,profiled);double total=ms(begin);counts=summary(mesh);
            if(counts.vertices!=expected[method].vertices || counts.quads!=expected[method].quads || counts.triangles!=expected[method].triangles)
                throw std::runtime_error("Profiling changed output counts");
            if(method==0 && dc.qefInteriorReturns+dc.qefBoundaryCells!=dc.activeCells)
                throw std::runtime_error("Inconsistent QEF counters");
            double profilePreparation=method==0?std::max(0.,total-dc.samplingMs-dc.hermiteMs-dc.qefMs-dc.connectivityMs):pv.preparationMs;
            // DC preparation is a residual including setup and unclassified overhead.
            out<<caseName<<','<<(method==0?"dc":"paraview-amr")<<','<<cellCount(volume.grid)<<','<<median(times[method])<<','<<*std::min_element(times[method].begin(),times[method].end())<<','<<*std::max_element(times[method].begin(),times[method].end())<<','<<total<<','<<profilePreparation<<','<<dc.samplingMs<<','<<dc.hermiteMs<<','<<dc.qefMs<<','<<dc.connectivityMs<<','<<pv.filterMs<<','<<pv.conversionMs<<','<<dc.activeCells<<','<<dc.edgeVisits<<','<<dc.uniqueEdges<<','<<dc.fieldQueries<<','<<dc.qefInteriorReturns<<','<<dc.qefBoundaryCells<<','<<dc.hermiteEvaluations<<','<<dc.hermiteCacheHits<<','<<counts.vertices<<','<<counts.quads<<','<<counts.triangles<<'\n';
            std::cout<<caseName<<' '<<(method==0?"dc":"paraview-amr")<<" median_ms="<<median(times[method])<<" vertices="<<counts.vertices<<'\n';
        }
        raw.flush();out.flush();if(!raw || !out)throw std::runtime_error("Report write failed");
    }
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
