#include "Field.h"
#include <cmath>
#include <limits>
#include <stdexcept>
namespace geometry {
std::size_t cellCount(const Grid& g) {
    std::size_t count=1;
    for(int n:g.cells) {
        if(n<=0 || n==std::numeric_limits<int>::max() || count>std::numeric_limits<std::size_t>::max()/static_cast<std::size_t>(n))
            throw std::invalid_argument("Invalid grid dimensions");
        count*=static_cast<std::size_t>(n);
    }
    return count;
}
void validateGrid(const Grid& g) {
    cellCount(g);
    if(!g.origin.allFinite() || !g.spacing.allFinite() || (g.spacing.array()<=0).any())
        throw std::invalid_argument("Invalid grid geometry");
    const Vec3 upper=g.origin+g.spacing.cwiseProduct(Vec3(g.cells[0],g.cells[1],g.cells[2]));
    if(!upper.allFinite() || (upper.array()<=g.origin.array()).any()) throw std::invalid_argument("Invalid grid bounds");
}
std::size_t cellIndex(const Grid& g,int x,int y,int z) {
    if(x<0 || y<0 || z<0 || x>=g.cells[0] || y>=g.cells[1] || z>=g.cells[2]) throw std::out_of_range("Cell index");
    return (static_cast<std::size_t>(x)*g.cells[1]+y)*g.cells[2]+z;
}
}
