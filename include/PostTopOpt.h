#pragma once
#include "SurfaceExtraction.h"
#include "ParaViewExtraction.h"
namespace geometry {
// One PostTopOpt API; each backend retains its own numerical settings.
enum class Backend { Dc, Emilk, Mc, CgalDc, ParaViewAmr };
struct ExtractionSettings {
    ContourOptions dc;
    CgalContourOptions cgal;
    ParaViewContourOptions paraview;
    double mcIsovalue=0;
};
Mesh extract(Backend,const Grid&,const Field&,const Gradient& = {},const ExtractionSettings& = {}); // метод построения поверхности, когда поле уже задано функцией
Mesh extract(Backend,const CellScalarVolume&,double densityThreshold,const ExtractionSettings& = {}); // метод построения поверхности, когда есть плотности в элементах
}
