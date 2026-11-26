#include "BuildIsoSurface.h"
#include "DualContouring.h"
#include <chrono>
#include <cmath>
#include <stdexcept>

namespace postopdc {
Surface buildIsoSurface(const std::vector<Vec3>& points,
    const std::vector<const Element*>& elements, const ExtractionOptions& options) {
    if (options.extractionSeconds) *options.extractionSeconds=0;
    if (!(options.iso>0 && options.iso<1) || !std::isfinite(options.regularization)
        || options.regularization<0) throw std::invalid_argument("Invalid DC options");
    // One shared VF per original node: arithmetic mean of incident element VF.
    std::vector<float> nodalVF(points.size(),0);
    std::vector<std::size_t> count(points.size(),0);
    for (const Element* element:elements) {
        if (!element || (element->type!=ElementType::Tetra4 &&
            element->type!=ElementType::Wedge6 && element->type!=ElementType::Hexa8))
            throw std::invalid_argument("Expected Tetra4, Wedge6 or Hexa8 element");
        if (!std::isfinite(element->vf) || element->vf<0 || element->vf>1)
            throw std::invalid_argument("VF must be in [0,1]");
        for (int a=0; a<element->nodeCount(); ++a) {
            const auto node=element->nodes[a];
            if (node>=points.size()) throw std::invalid_argument("Invalid node index");
            nodalVF[node]+=element->vf;
            ++count[node];
        }
    }
    for (std::size_t i=0; i<points.size(); ++i)
        if (count[i]) nodalVF[i]/=static_cast<float>(count[i]);
    // Original mesh only: no transfer, regular grid, padding or numerical gradient.
    const auto begin=options.extractionSeconds ? std::chrono::steady_clock::now()
                                             : std::chrono::steady_clock::time_point{};
    auto surface=dualContour(points,elements,nodalVF,options);
    if (options.extractionSeconds)
        *options.extractionSeconds=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-begin).count();
    return surface;
}
} // namespace postopdc
