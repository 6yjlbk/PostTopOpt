#include "DualContouring.h"
#include "ObjWriter.h"
#include <iostream>
int main(int argc,char** argv) {
    try {
        geometry::Grid grid;
        auto mesh=geometry::dualContour(grid,[](const geometry::Vec3& p){ return p.norm()-2; },
            [](const geometry::Vec3& p)->geometry::Vec3 { return p.normalized(); });
        const auto path=argc>1 ? argv[1] : "sphere.obj";
        geometry::writeObj(mesh,path);
        std::cout<<mesh.vertices.size()<<" vertices, "<<mesh.quads.size()<<" quads -> "<<path<<'\n';
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
