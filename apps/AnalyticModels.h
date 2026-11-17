#pragma once
#include "CellScalarVolume.h"
#include <string>
namespace application {
geometry::CellScalarVolume generateAnalytic(const std::string&, int cells, int complexity);
}
