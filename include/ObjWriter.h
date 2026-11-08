#pragma once
#include "Mesh.h"
#include <filesystem>
namespace geometry {
void writeObj(const Mesh&, const std::filesystem::path&);
}
