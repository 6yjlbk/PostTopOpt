#include "ExampleRunner.h"
#include <algorithm>

int main(int argc, char** argv) {
    // A hard union: intersection seams remain sharp, without smooth blending.
    const auto field = [](const geometry::Vec3& p) {
        return std::min({
            (p - geometry::Vec3(-0.55, 0, 0)).norm() - 0.9,
            (p - geometry::Vec3( 0.55, 0, 0)).norm() - 0.9,
            (p - geometry::Vec3(0, 0.65, 0.25)).norm() - 0.75
        });
    };
    return examples::run(argc, argv, "sphere_union.obj", examples::grid(1.8, 60), field);
}

