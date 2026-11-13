#include "AnalyticModels.h"
#include "VolumeIO.h"
#include "PostTopOpt.h"
#include "ObjWriter.h"
#ifdef POSTTOPOPT_APP_WITH_VTK
#include "PlyDistanceField.h"
#endif
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <locale>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace geometry;
using Clock=std::chrono::steady_clock;
double elapsed(Clock::time_point start) {
    return std::chrono::duration<double,std::milli>(Clock::now()-start).count();
}
bool boolean(const std::string& text) {
    if(text=="true" || text=="1") return true;
    if(text=="false" || text=="0") return false;
    throw std::invalid_argument("Boolean requires true|false or 1|0");
}
int integer(const std::string& text) {
    std::size_t used=0;int result=std::stoi(text,&used);
    if(used!=text.size())throw std::invalid_argument("Invalid integer");return result;
}
double number(const std::string& text) {
    std::size_t used=0;double result=std::stod(text,&used);
    if(used!=text.size() || !std::isfinite(result))throw std::invalid_argument("Invalid number");return result;
}
struct Method {std::string name;Backend backend;};
Method method(const std::string& name) {
    if(name=="dc")return {name,Backend::Dc};
    if(name=="emilk")return {name,Backend::Emilk};
    if(name=="mc")return {name,Backend::Mc};
    if(name=="cgal-dc")return {name,Backend::CgalDc};
    if(name=="paraview-amr")return {name,Backend::ParaViewAmr};
    throw std::invalid_argument("Unknown method: "+name);
}
struct Counts {std::size_t vertices,quads,triangles;bool operator==(const Counts&) const = default;};
Counts check(const Mesh& mesh) {
    if(mesh.vertices.empty() || (mesh.quads.empty() && mesh.triangles.empty()))throw std::runtime_error("Empty mesh");
    for(const auto& p:mesh.vertices)if(!p.allFinite())throw std::runtime_error("Non-finite vertex");
    for(const auto& face:mesh.quads)for(auto i:face)if(i>=mesh.vertices.size())throw std::runtime_error("Invalid quad index");
    for(const auto& face:mesh.triangles)for(auto i:face)if(i>=mesh.vertices.size())throw std::runtime_error("Invalid triangle index");
    return {mesh.vertices.size(),mesh.quads.size(),mesh.triangles.size()};
}
double median(std::vector<double> values) {
    std::sort(values.begin(),values.end());return (values[(values.size()-1)/2]+values[values.size()/2])/2.;
}
}
int main(int argc,char** argv) {try {
    const auto processStart=Clock::now();
    std::string input,selected="dc";
    std::filesystem::path output,saveVolume;
    int cells=64,complexity=4,repeats=5;
    bool benchmark=true,exportObj=true,cellsExplicit=false,complexityExplicit=false;
    double widthCells=2.,threshold=.5;
    ExtractionSettings settings;
    for(int i=1;i<argc;++i) {
        const std::string arg=argv[i];
        if(arg=="--help") {
            std::cout<<"PostTopOptRun --input file.ply|file.csv|file.ptov|model --output-dir directory\n"
                "Models: sphere torus sphere-union drilled-box cone gyroid channels lattice shells\n"
                "--cells 64 --complexity 4 --width-cells 2 --threshold 0.5\n"
                "--method dc|emilk|mc|cgal-dc|paraview-amr|all|compare\n"
                "--benchmark true|false (default: true) --repeats 5 --export true|false --save-volume path.csv\n"
                "--lambda 0.001 --regularization true --constrained true --bisection true\n"
                "--cache true --qef true\n";return 0;
        }
        if(++i==argc)throw std::invalid_argument("Missing option value: "+arg);
        const std::string value=argv[i];
        if(arg=="--input")input=value;
        else if(arg=="--output-dir")output=value;
        else if(arg=="--save-volume")saveVolume=value;
        else if(arg=="--method")selected=value;
        else if(arg=="--cells"){cells=integer(value);cellsExplicit=true;}
        else if(arg=="--complexity"){complexity=integer(value);complexityExplicit=true;}
        else if(arg=="--repeats")repeats=integer(value);
        else if(arg=="--width-cells")widthCells=number(value);
        else if(arg=="--threshold")threshold=number(value);
        else if(arg=="--benchmark")benchmark=boolean(value);
        else if(arg=="--export")exportObj=boolean(value);
        else if(arg=="--lambda")settings.dc.regularization=number(value);
        else if(arg=="--regularization")settings.dc.regularizationEnabled=boolean(value);
        else if(arg=="--constrained")settings.dc.constrainToCell=settings.cgal.constrainToCell=boolean(value);
        else if(arg=="--bisection")settings.dc.bisection=settings.cgal.bisection=boolean(value);
        else if(arg=="--cache")settings.dc.cacheHermite=boolean(value);
        else if(arg=="--qef")settings.dc.optimizeQef=boolean(value);
        else throw std::invalid_argument("Unknown option: "+arg);
    }
    if(input.empty() || output.empty())throw std::invalid_argument("--input and --output-dir required");
    if(cells<16 || cells>500 || complexity<1 || complexity>12 || (benchmark && (repeats<2 || repeats>50)) || widthCells<=0 || widthCells>4 || threshold<=0 || threshold>=1 || settings.dc.regularization<0)
        throw std::invalid_argument("Expected cells:16..500 complexity:1..12 repeats:2..50 width-cells:(0,4] threshold:(0,1) lambda>=0");
    std::vector<Method> methods;
    if(selected=="all")for(const char* name:{"dc","emilk","mc","cgal-dc","paraview-amr"})methods.push_back(method(name));
    else if(selected=="compare")methods={method("dc"),method("paraview-amr")};
    else methods.push_back(method(selected));
    for(const auto& m:methods) {
        if(m.backend==Backend::CgalDc && !cgalAvailable())throw std::runtime_error("CGAL backend disabled");
        if(m.backend==Backend::ParaViewAmr && !paraViewAvailable())throw std::runtime_error("ParaView backend disabled");
    }
    if(std::filesystem::exists(output) && (!std::filesystem::is_directory(output) || !std::filesystem::is_empty(output)))
        throw std::runtime_error("Output directory must be new or empty");
    if(!saveVolume.empty() && (saveVolume.extension()!=".csv" || std::filesystem::exists(saveVolume)))
        throw std::runtime_error("--save-volume requires a new .csv path");
    CellScalarVolume volume;
    const auto extension=std::filesystem::path(input).extension().string();
    std::string sourceKind;
    if(extension==".ply") {
        if(complexityExplicit)throw std::invalid_argument("--complexity applies only to analytic models");
#ifdef POSTTOPOPT_APP_WITH_VTK
        sourceKind="ply-synthetic-volume-fraction";
        auto model=makePlyDistanceField(input,cells);
        const double width=model.grid.spacing.minCoeff()*widthCells;
        volume=sampleCells(model.grid,[&](const Vec3& p){return std::clamp(.5-model.field(p)/width,0.,1.);});
#else
        throw std::runtime_error("PLY requires POSTTOPOPT_WITH_VTK=ON");
#endif
    } else if(extension==".csv" || extension==".ptov") {
        if(cellsExplicit || complexityExplicit)throw std::invalid_argument("Saved volume already defines its grid; omit --cells and --complexity");
        sourceKind="saved-cell-volume-fraction";volume=application::loadVolume(input);
    } else if(extension.empty()) {
        sourceKind="analytic-synthetic-volume-fraction";volume=application::generateAnalytic(input,cells,complexity);
    } else throw std::invalid_argument("Input must be a model name, .ply, .csv or .ptov");
    volume.validate();
    const double preparationVfMs=elapsed(processStart);
    settings.dc.gradientStep=settings.cgal.gradientStep=volume.grid.spacing.minCoeff()*1e-3;
    std::filesystem::create_directories(output);
    if(!saveVolume.empty())application::saveVolume(volume,saveVolume);
    std::ofstream config(output/"settings.txt"),raw,report;
    if(benchmark) {
        raw.open(output/"runs.csv");report.open(output/"summary.csv");
        if(!raw || !report)throw std::runtime_error("Cannot create benchmark reports");
    }
    config<<std::setprecision(17)<<"input="<<input<<"\nsource="<<sourceKind<<"\nmethod="<<selected
        <<"\nbenchmark="<<benchmark<<"\nrepeats="<<(benchmark?repeats:1)<<"\nexport="<<exportObj
        <<"\ncells="<<volume.grid.cells[0]<<','<<volume.grid.cells[1]<<','<<volume.grid.cells[2]
        <<"\norigin="<<volume.grid.origin.transpose()<<"\nspacing="<<volume.grid.spacing.transpose()
        <<"\ncomplexity="<<complexity<<"\nwidth_cells="<<widthCells<<"\nthreshold="<<threshold
        <<"\nlambda="<<settings.dc.regularization<<"\nregularization="<<settings.dc.regularizationEnabled
        <<"\nconstrained="<<settings.dc.constrainToCell<<"\nbisection="<<settings.dc.bisection
        <<"\ngradient_step="<<settings.dc.gradientStep<<"\ncache="<<settings.dc.cacheHermite<<"\nqef="<<settings.dc.optimizeQef
        <<"\nsave_volume="<<saveVolume.string()<<"\n";
    if(!config)throw std::runtime_error("Cannot create settings");
    if(benchmark) {
        raw.imbue(std::locale::classic());report.imbue(std::locale::classic());
        raw<<std::fixed<<std::setprecision(3)<<"method,repeat,isosurf_s,vertices,quads,triangles\n";
        // UTF-8 BOM lets Excel detect the Russian metric descriptions.
        report<<"\xEF\xBB\xBF"<<"metric,description";
        for(const auto& m:methods)report<<','<<m.name;
        report<<'\n';
    }
    struct ReportRow {std::string key,label;std::vector<std::string> values;};
    std::vector<ReportRow> rows;
    auto record=[&](std::size_t methodIndex,const std::string& key,const std::string& label,const std::string& value) {
        auto it=std::find_if(rows.begin(),rows.end(),[&](const auto& row){return row.key==key;});
        if(it==rows.end()){rows.push_back({key,label,std::vector<std::string>(methods.size())});it=std::prev(rows.end());}
        it->values[methodIndex]=value;
    };
    auto seconds=[](double milliseconds) {
        std::ostringstream text;text.imbue(std::locale::classic());
        text<<std::fixed<<std::setprecision(3)<<milliseconds/1000.;return text.str();
    };
    std::vector<std::vector<double>> times(methods.size());
    std::vector<Counts> expected(methods.size());
    std::vector<double> exportTimes(methods.size(),0.);
    if(benchmark)for(const auto& m:methods)check(extract(m.backend,volume,threshold,settings));
    // Rotate method order between repetitions; mesh destruction and checking are outside timing.
    for(int rep=0;rep<(benchmark?repeats:1);++rep)for(std::size_t order=0;order<methods.size();++order) {
        const auto i=(order+rep)%methods.size();
        const auto start=Clock::now();const auto mesh=extract(methods[i].backend,volume,threshold,settings);
        const double extractionMs=elapsed(start);const auto counts=check(mesh);
        if(rep && !(counts==expected[i]))throw std::runtime_error("Output counts changed between repetitions");
        expected[i]=counts;times[i].push_back(extractionMs);
        if(!benchmark && exportObj) {
            const auto exportStart=Clock::now();writeObj(mesh,output/(methods[i].name+".obj"));exportTimes[i]=elapsed(exportStart);
        }
        if(benchmark)raw<<methods[i].name<<','<<rep<<','<<extractionMs/1000.<<','<<counts.vertices<<','<<counts.quads<<','<<counts.triangles<<'\n';
    }
    // Diagnostic timings are measured separately; normal repetitions have no profiling overhead.
    for(std::size_t i=0;i<methods.size();++i) {
        DcProfile dc;ParaViewProfile pv;auto diagnostic=settings;
        diagnostic.dc.profile=&dc;diagnostic.paraview.profile=&pv;
        double exportMs=exportTimes[i],profileMs=0;
        if(benchmark) {
            const auto start=Clock::now();
            const auto mesh=extract(methods[i].backend,volume,threshold,diagnostic);
            profileMs=elapsed(start);
            if(!(check(mesh)==expected[i]))throw std::runtime_error("Diagnostic/export output counts changed");
            if(exportObj){const auto exportStart=Clock::now();writeObj(mesh,output/(methods[i].name+".obj"));exportMs=elapsed(exportStart);}
        }
        const auto counts=expected[i];const double med=median(times[i]);
        if(benchmark) {
            auto time=[&](const char* key,const char* label,double value){record(i,key,label,seconds(value));};
            time("preparation_vf_s","Время подготовки VF",preparationVfMs);
            time("isosurf_median_s","Медианное время построения изоповерхности",med);
            record(i,"isosurf_min_s - isosurf_max_s","Минимальное — максимальное время построения изоповерхности",
                seconds(*std::min_element(times[i].begin(),times[i].end()))+" - "+seconds(*std::max_element(times[i].begin(),times[i].end())));
            time("build_isosurf_s","Время подготовки VF и построения изоповерхности",preparationVfMs+med);
            time("algorithm_without_vf_s","Общее время алгоритма без подготовки VF",med);
            time("algorithm_with_vf_s","Общее время алгоритма с подготовкой VF",preparationVfMs+med);
            record(i,"export_s","Время экспорта OBJ",exportObj?seconds(exportMs):"");
            record(i,"vertices","Количество вершин",std::to_string(counts.vertices));
            record(i,"quads","Количество четырёхугольников",std::to_string(counts.quads));
            record(i,"triangles","Количество треугольников",std::to_string(counts.triangles));
            time("profile_isosurf_s","Время построения изоповерхности с профилированием",profileMs);
            const bool ownDc=methods[i].backend==Backend::Dc || methods[i].backend==Backend::Emilk;
            const bool vtk=methods[i].backend==Backend::ParaViewAmr;
            auto stage=[&](const char* key,const char* label,bool applicable,double value){record(i,key,label,applicable?seconds(value):"");};
            stage("sampling_s","Время вычисления VF в узлах",ownDc,dc.samplingMs);
            stage("hermite_s","Время подготовки данных Эрмита",ownDc,dc.hermiteMs);
            stage("qef_s","Время решения QEF",ownDc,dc.qefMs);
            stage("connectivity_s","Время построения связности",ownDc,dc.connectivityMs);
            stage("vtk_preparation_s","Время подготовки данных VTK",vtk,pv.preparationMs);
            stage("filter_s","Время работы фильтра ParaView",vtk,pv.filterMs);
            stage("conversion_s","Время преобразования результата ParaView в Mesh",vtk,pv.conversionMs);
            auto count=[&](const char* key,const char* label,std::size_t value){record(i,key,label,ownDc?std::to_string(value):"");};
            count("field_queries","Количество запросов поля VF",dc.fieldQueries);
            count("active_cells","Количество активных ячеек",dc.activeCells);
            count("edge_visits","Количество посещений пересечённых рёбер",dc.edgeVisits);
            count("unique_edges","Количество уникальных пересечённых внутренних рёбер",dc.uniqueEdges);
            count("hermite_evaluations","Количество вычислений данных Эрмита",dc.hermiteEvaluations);
            count("hermite_cache_hits","Количество попаданий в кэш Эрмита",dc.hermiteCacheHits);
            count("qef_interior","Количество ранних возвратов QEF",dc.qefInteriorReturns);
            count("qef_boundary","Количество граничных поисков QEF",dc.qefBoundaryCells);
        }
        std::cout<<methods[i].name<<" isosurf_s="<<seconds(med)<<" vertices="<<counts.vertices<<'\n';
    }
    if(benchmark)for(const auto& row:rows) {
        report<<row.key<<','<<std::quoted(row.label);
        for(const auto& value:row.values)report<<','<<value;
        report<<'\n';
    }
    config.flush();if(!config)throw std::runtime_error("Settings write failed");
    if(benchmark) {
        raw.flush();report.flush();
        if(!raw || !report)throw std::runtime_error("Benchmark report write failed");
    }
    std::cout<<"preparation_vf_s="<<seconds(preparationVfMs)<<" output="<<output.string()<<'\n';
} catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}}
