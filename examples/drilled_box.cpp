#include "ExampleRunner.h"
#include <algorithm>
#include <cmath>

int main(int argc, char** argv) {
    // Subtract three infinite, perpendicular cylinders from a box.
    const auto field = [](const geometry::Vec3& p) {
        const double box = p.cwiseAbs().maxCoeff() - 1.15;
        const double channels = std::min({
            std::hypot(p.y(), p.z()) - 0.46,
            std::hypot(p.x(), p.z()) - 0.46,
            std::hypot(p.x(), p.y()) - 0.46
        });
        return std::max(box, -channels);
    };
    return examples::run(argc, argv, "drilled_box.obj", examples::grid(1.5, 60), field);
}

