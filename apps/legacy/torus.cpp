#include "ExampleRunner.h"
#include <cmath>

int main(int argc, char** argv) {
    // Major radius: 1.05; tube radius: 0.38; symmetry axis: z.
    const auto field = [](const geometry::Vec3& p) {
        return std::hypot(std::hypot(p.x(), p.y()) - 1.05, p.z()) - 0.38;
    };
    return examples::run(argc, argv, "torus.obj", examples::grid(1.8, 60), field);
}

