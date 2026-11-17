#include <utility>
#include "CellScalarVolume.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace geometry {
void CellScalarVolume::validate() const {
    validateGrid(grid);
    if(values.size()!=cellCount(grid)) throw std::invalid_argument("Cell value count differs from grid");
    for(double v:values) if(!std::isfinite(v)) throw std::invalid_argument("Non-finite cell value");
}
Field linearCellField(std::shared_ptr<const CellScalarVolume> volume) {
    if(!volume) throw std::invalid_argument("Null cell volume");
    volume->validate();
    return [volume](const Vec3& p) {
        if(!p.allFinite()) throw std::invalid_argument("Non-finite query point");
        const auto& g=volume->grid;
        int lo[3],hi[3]; double t[3];
        for(int d=0;d<3;++d) {
            const double u=std::clamp((p[d]-g.origin[d])/g.spacing[d]-0.5,0.,static_cast<double>(g.cells[d]-1));
            lo[d]=static_cast<int>(std::floor(u)); hi[d]=std::min(lo[d]+1,g.cells[d]-1); t[d]=u-lo[d];
        }
        double value=0;
        for(int mask=0;mask<8;++mask) {
            int at[3]; double weight=1;
            for(int d=0;d<3;++d) { bool upper=(mask&(1<<d))!=0; at[d]=upper?hi[d]:lo[d]; weight*=upper?t[d]:1-t[d]; }
            value+=weight*volume->values[(static_cast<std::size_t>(at[0])*g.cells[1]+at[1])*g.cells[2]+at[2]];
        }
        return value;
    };
}
Field densityField(std::shared_ptr<const CellScalarVolume> volume,double threshold) {
    if(!std::isfinite(threshold)) throw std::invalid_argument("Non-finite density threshold");
    auto rho=linearCellField(std::move(volume));
    return [rho,threshold](const Vec3& p){return threshold-rho(p);};
}
CellScalarVolume sampleCells(const Grid& g,const Field& field) {
    validateGrid(g); if(!field) throw std::invalid_argument("Missing field");
    CellScalarVolume result{g,std::vector<double>(cellCount(g))};
    for(int x=0;x<g.cells[0];++x) for(int y=0;y<g.cells[1];++y) for(int z=0;z<g.cells[2];++z) {
        Vec3 p=g.origin+g.spacing.cwiseProduct(Vec3(x+0.5,y+0.5,z+0.5));
        result.values[cellIndex(g,x,y,z)]=field(p);
    }
    result.validate(); return result;
}
}
