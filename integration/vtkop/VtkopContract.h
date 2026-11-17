#pragma once
// PROPOSED CONTRACT ONLY, NOT THE EXISTING vtkop API.
// Implement/map these operations inside the real vtkop on the target machine.
// This file deliberately has no VTK includes and no backend implementation.
// Replace this include boundary with the real vtkop header during integration.
#include <array>
#include <cstddef>
#include <memory>
#include <vector>

namespace postop::vtkop_contract {
struct UniformGridObject;
struct AmrObject;
struct ContourObject;
struct PolyDataObject;
using UniformGrid = std::shared_ptr<UniformGridObject>;
using Amr = std::shared_ptr<AmrObject>;
using Contour = std::shared_ptr<ContourObject>;
using PolyData = std::shared_ptr<PolyDataObject>;

// Handles retain their underlying VTK objects; destruction must not throw.
// Operations throw std::runtime_error on VTK failure, never return null
// for successful creation. Empty geometry is represented by a valid object.
// Validate all size_t/int conversions to vtkIdType in the implementation.
UniformGrid createUniformGrid(
    const std::array<double, 3>& origin,
    const std::array<double, 3>& spacing,
    const std::array<int, 3>& pointDimensions);

// Copy before returning, set active cell scalars; do not retain values pointer.
// x-fastest order, double precision, no scaling of VF to [0,255].
void setCellScalars(const UniformGrid& grid, const char* name,
                    const double* values, std::size_t count);

// One level, one block at index 0; retain block ownership.
Amr createSingleBlockAmr(const UniformGrid& grid);
Contour createAmrDualContour();
void setInput(const Contour& contour, const Amr& input);
void setCellArray(const Contour& contour, const char* name);
void setIsovalue(const Contour& contour, double value);
void setCapping(const Contour& contour, bool enabled);
void setMergePoints(const Contour& contour, bool enabled);

// Use the serial controller behavior of src/ParaViewExtraction.cpp.
// Restore the global controller on success AND exceptions; coordinate with
// the host's other controller users. A private mutex alone is insufficient.
// Attach VTK ErrorEvent observation and report failed Update calls.
void update(const Contour& contour);

// Traverse nonempty vtkPolyData leaves in their original order, skipping
// other leaf types as the existing implementation does. Returned handles
// remain valid after the producer handle is released.
std::vector<PolyData> polyDataBlocks(const Contour& contour);

// vtkTriangleFilter, PassLinesOff/PassVertsOff; retain resulting data.
// No smoothing, decimation, welding or orientation changes.
PolyData triangulate(const PolyData& input);
std::size_t pointCount(const PolyData& data);
std::size_t triangleCount(const PolyData& data);
std::array<double, 3> point(const PolyData& data, std::size_t index);
// Reject non-triangles and negative VTK IDs before conversion to size_t.
// Preserve original face order and orientation. Access should not rescan
// the entire cell array for every index (use VTK indexed access or a cache).
std::array<std::size_t, 3> triangle(const PolyData& data, std::size_t index);
} // namespace postop::vtkop_contract
