#include "PostTopOpt.h"
#include <memory>
#include <stdexcept>
namespace geometry {
Mesh extract(Backend backend,const Grid& grid,const Field& field,const Gradient& gradient,const ExtractionSettings& settings) {
    switch(backend) {
    case Backend::Dc:return dualContourGeneral(grid,field,gradient,settings.dc);
    case Backend::Emilk:return dualContourEmilk(grid,field,gradient,settings.dc);
    case Backend::Mc:return marchingCubes(grid,field,settings.mcIsovalue);
    case Backend::CgalDc:return cgalDualContour(grid,field,gradient,settings.cgal);
    case Backend::ParaViewAmr:throw std::invalid_argument("ParaView AMR requires cell volume fractions; sample an explicit density field with sampleCells");
    }
    throw std::invalid_argument("Unknown backend");
}
Mesh extract(Backend backend,const CellScalarVolume& volume,double iso,const ExtractionSettings& settings) {

    if(backend==Backend::ParaViewAmr) {
        auto options=settings.paraview; options.isovalue=iso;
        return paraViewDualContour(volume,options);
    }
    // Synchronous extraction never retains Field: borrow the immutable input instead of copying it.
    auto field=densityField(std::shared_ptr<const CellScalarVolume>(&volume,[](const CellScalarVolume*){}),iso);
    auto options=settings;
    // densityField already incorporates the physical density iso.
    options.dc.isovalue=0; options.cgal.isovalue=0; options.mcIsovalue=0;
    return extract(backend,volume.grid,field,{},options);
}
}

