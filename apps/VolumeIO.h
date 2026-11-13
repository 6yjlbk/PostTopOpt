#pragma once
#include "CellScalarVolume.h"
#include <filesystem>
namespace application {
geometry::CellScalarVolume loadVolume(const std::filesystem::path&);
void saveVolume(const geometry::CellScalarVolume&, const std::filesystem::path&);
}
