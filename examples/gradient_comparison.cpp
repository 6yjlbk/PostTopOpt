#include "ExampleRunner.h"
#include "CgalExtraction.h"
#include <algorithm>
#include <cmath>

// Compare only the gradient source; both modes use CGAL DC and биссекция.
int main(int argc, char** argv) {
    try {
        if (argc != 3) throw std::invalid_argument("Usage: PostTopOptGradientComparison sphere|cone numeric|analytic");
        const std::string shape=argv[1], mode=argv[2];
        if ((shape!="sphere" && shape!="cone") || (mode!="numeric" && mode!="analytic"))
            throw std::invalid_argument("Expected sphere|cone and numeric|analytic");
        geometry::Field field;
        geometry::Gradient analytic;
        auto grid=examples::grid(3.0,80);
        if (shape=="sphere") {
            field=[](const geometry::Vec3& p) {return p.norm()-2.0;};
            analytic=[](const geometry::Vec3& p)->geometry::Vec3 {
                const double r=p.norm();
                if (r==0) throw std::invalid_argument("Sphere gradient undefined at center");
                return p/r; // gradient
            };
        } else { //Cone
            // Solid cone: R=1, H=2, base z=0, apex z=2.
            // f=max(r + (R/H)*z - R, -z). f<0 implies 0<z<H,
            // so a separate z-H cap is redundant and is deliberately omitted.
            // This field is not a signed distance function.
            grid=examples::grid(1.5,80);
            grid.origin.z()=-0.5;
            field=[](const geometry::Vec3& p) {
                return std::max(std::hypot(p.x(),p.y())+0.5*p.z()-1.0,-p.z());
            };
            analytic=[](const geometry::Vec3& p)->geometry::Vec3 {
                const double r=std::hypot(p.x(),p.y());
                const double side=r+0.5*p.z()-1.0;
                // At the base rim the derivative is not unique: choose the base.
                if (-p.z()>=side) return geometry::Vec3(0,0,-1);
                // On the axis the radial derivative is undefined. Choose the
                // symmetric subgradient, consistent with central differences.
                if (r==0) return geometry::Vec3(0,0,0.5);
                return geometry::Vec3(p.x()/r,p.y()/r,0.5);
            };
        }
        geometry::CgalContourOptions options;
        options.bisection=false; // Биссекция in both comparisons.
        options.isovalue=0;
        options.constrainToCell=false;
        options.gradientStep=1e-5;
        const geometry::Gradient gradient=mode=="analytic" ? analytic : geometry::Gradient{};
        const auto begin=std::chrono::steady_clock::now();
        const auto mesh=geometry::cgalDualContour(grid,field,gradient,options);
        const auto elapsed=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-begin).count();
        const auto path=std::filesystem::path(POSTTOPOPT_RESULTS_DIR)/(shape+"_cgal-dc_"+mode+".obj");
        std::filesystem::create_directories(path.parent_path());
        geometry::writeObj(mesh,path);
        std::cout<<shape<<" "<<mode<<": "<<mesh.vertices.size()<<" vertices, "<<mesh.quads.size()
                 <<" quads, extraction_ms="<<elapsed<<" -> "<<path.string()<<'\n';
        return 0;
    } catch (const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
