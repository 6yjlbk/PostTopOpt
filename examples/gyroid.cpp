#include "ExampleRunner.h"
#include <algorithm>
#include <cmath>
#include <numbers>

int main(int argc, char** argv) {
    // A finite sheet: |g| < 0.28, intersected with a box.
    // This is an implicit band, not a constant physical wall thickness.
    const auto field = [](const geometry::Vec3& p) {
        const geometry::Vec3 q = std::numbers::pi * p;
        const double g = std::sin(q.x()) * std::cos(q.y())
                       + std::sin(q.y()) * std::cos(q.z())
                       + std::sin(q.z()) * std::cos(q.x());
        const double sheet = std::abs(g) - 0.28;
        const double box = p.cwiseAbs().maxCoeff() - 1.35;
        return std::max(sheet, box);
    };
    return examples::run(argc, argv, "gyroid.obj", examples::grid(1.6, 80), field);
}

