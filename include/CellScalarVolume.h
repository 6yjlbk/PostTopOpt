#pragma once
#include "Field.h"
#include <memory>
namespace geometry {
// Values are CELL-associated, not nodal. Index: (x*Ny+y)*Nz+z.
struct CellScalarVolume {
    Grid grid;
    std::vector<double> values;
    void validate() const;
};
// Trilinear interpolation on cell centers; constant continuation outside their hull.
// This does not fill/cap the physical domain and is not a shape-function interpolation.
Field linearCellField(std::shared_ptr<const CellScalarVolume>); //трилинейная интерполяция виртуальной плотности между центрами КЭ
// Converts high density to negative-inside convention, with zero surface threshold.
Field densityField(std::shared_ptr<const CellScalarVolume>, double threshold); //учет порога материала
CellScalarVolume sampleCells(const Grid&, const Field&);
}
