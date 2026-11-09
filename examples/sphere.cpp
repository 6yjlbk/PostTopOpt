#include "ExampleRunner.h"
int main(int argc, char** argv) {
    geometry::Grid grid;
    return examples::run(argc, argv, "sphere.obj", grid,
        [](const geometry::Vec3& p) { return p.norm()-2; },
        [](const geometry::Vec3& p)->geometry::Vec3 { return p.normalized(); });
}
