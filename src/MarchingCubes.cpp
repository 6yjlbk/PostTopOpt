#include "MarchingCubes.h"
#include "MarchingCubesTable.h"
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace geometry {
namespace {
constexpr int corners[8][3] = {
    {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
    {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}
};
constexpr int edges[12][2] = {
    {0,1}, {1,2}, {2,3}, {3,0}, {4,5}, {5,6},
    {6,7}, {7,4}, {0,4}, {1,5}, {2,6}, {3,7}
};
}
Mesh marchingCubes(const Grid& grid, const Field& field, double isovalue) {
    if (!field || !grid.origin.allFinite() || !grid.spacing.allFinite() ||
        (grid.spacing.array() <= 0).any() || !std::isfinite(isovalue))
        throw std::invalid_argument("Invalid field, grid or isovalue");
    auto nodes = grid.cells;
    std::size_t nodeCount = 1;
    for (auto& count : nodes) {
        if (count <= 0 || count == std::numeric_limits<int>::max())
            throw std::invalid_argument("Invalid grid dimensions");
        ++count;
        if (nodeCount > std::numeric_limits<std::size_t>::max() / static_cast<std::size_t>(count))
            throw std::invalid_argument("Grid dimensions overflow");
        nodeCount *= static_cast<std::size_t>(count);
    }
    const Vec3 maximum = grid.origin + grid.spacing.cwiseProduct(
        Vec3(grid.cells[0], grid.cells[1], grid.cells[2]));
    if (!maximum.allFinite()) throw std::invalid_argument("Grid coordinates overflow");
    const auto index = [&](int x, int y, int z) {
        return (static_cast<std::size_t>(x) * nodes[1] + y) * nodes[2] + z;
    };
    const auto position = [&](int x, int y, int z) -> Vec3 {
        return grid.origin + grid.spacing.cwiseProduct(Vec3(x, y, z));
    };
    std::vector<double> values(nodeCount);
    for (int x=0; x<nodes[0]; ++x) for (int y=0; y<nodes[1]; ++y) for (int z=0; z<nodes[2]; ++z) {
        const double value = field(position(x,y,z)) - isovalue;
        if (!std::isfinite(value)) throw std::invalid_argument("Non-finite field value");
        values[index(x,y,z)] = value;
    }
    Mesh mesh;
    // Global grid-edge keys share vertices between cells. (node,node) also
    // welds intersections exactly at a grid node across different edges.
    std::map<std::pair<std::size_t,std::size_t>,std::size_t> vertexIds;
    for (int x=0; x<grid.cells[0]; ++x) for (int y=0; y<grid.cells[1]; ++y)
    for (int z=0; z<grid.cells[2]; ++z) {
        std::array<std::size_t,8> ids;
        std::array<Vec3,8> positions;
        std::array<double,8> f;
        int mask = 0;
        for (int c=0; c<8; ++c) {
            const int cx=x+corners[c][0], cy=y+corners[c][1], cz=z+corners[c][2];
            ids[c]=index(cx,cy,cz); positions[c]=position(cx,cy,cz); f[c]=values[ids[c]];
            if (f[c]<0) mask |= 1<<c;
        }
        const auto& row = detail::mcTriangles[mask];
        if (row[0]<0) continue;
        std::array<std::size_t,12> localIds{};
        std::array<bool,12> ready{};
        const auto vertex = [&](int edge) -> std::size_t {
            if (ready[edge]) return localIds[edge];
            int a=edges[edge][0], b=edges[edge][1];
            if (ids[a]>ids[b]) std::swap(a,b); // Deterministic interpolation order.
            std::pair<std::size_t,std::size_t> key{ids[a],ids[b]};
            if (f[a]==0) key={ids[a],ids[a]};
            else if (f[b]==0) key={ids[b],ids[b]};
            auto found=vertexIds.find(key);
            std::size_t id;
            if (found!=vertexIds.end()) id=found->second;
            else {
                const double scale=std::max(std::abs(f[a]),std::abs(f[b]));
                const double t=(f[a]/scale)/((f[a]/scale)-(f[b]/scale));
                const Vec3 p=positions[a]+t*(positions[b]-positions[a]);
                if (!p.allFinite()) throw std::runtime_error("Invalid MC intersection");
                id=mesh.vertices.size();
                mesh.vertices.push_back(p);
                vertexIds.emplace(key,id);
            }
            ready[edge]=true; localIds[edge]=id;
            return id;
        };
        for (int t=0; t<16 && row[t]>=0; t+=3) {
            // Table winds towards negative values; reverse for outward faces.
            std::array<std::size_t,3> triangle{vertex(row[t]),vertex(row[t+2]),vertex(row[t+1])};
            if (triangle[0]==triangle[1] || triangle[1]==triangle[2] || triangle[2]==triangle[0]) continue;
            const Vec3 normal=(mesh.vertices[triangle[1]]-mesh.vertices[triangle[0]])
                .cross(mesh.vertices[triangle[2]]-mesh.vertices[triangle[0]]);
            if (normal.squaredNorm()==0) continue;
            mesh.triangles.push_back(triangle);
        }
    }
    // Remove vertices used only by collapsed triangles at exact isovalues.
    const auto missing=std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> remap(mesh.vertices.size(),missing);
    for (const auto& t:mesh.triangles) for (auto id:t) remap[id]=0;
    std::vector<Vec3> compact;
    for (std::size_t i=0;i<remap.size();++i) if (remap[i]!=missing) {
        remap[i]=compact.size(); compact.push_back(mesh.vertices[i]);
    }
    for (auto& t:mesh.triangles) for (auto& id:t) id=remap[id];
    mesh.vertices=std::move(compact);
    return mesh;
}
}
