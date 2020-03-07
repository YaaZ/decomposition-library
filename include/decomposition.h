#pragma once


#include "decomposition-types.h"


namespace decomposition {


    /**
     * Finds all intersections of polygons (including self-intersections)
     * @param vertices vertex coordinates, new vertices will be appended here
     * @param polygonVertexIndices vertex indices per polygon
     * @return vertex linkage objects, one per resulting vertex. All linkage objects constitutes one or more Eulerian graphs
     */
    std::vector<VertexLinkage> insertSteinerVerticesForPolygons(std::vector<glm::dvec2>& vertices,
                                                                const std::vector<std::vector<int>>& polygonVertexIndices);

    /**
     * Transforms polygon vertex graph into set of non-intersecting polygons, keeping winding order of original polygons
     * @param vertices vertex coordinates
     * @param verticesLinkage vertex linkage objects, consumed by this function, see insertSteinerVerticesForPolygons()
     * @return generated non-intersecting polygons
     */
    std::vector<Polygon> decomposePolygonGraph(const std::vector<glm::dvec2>& vertices,
                                               std::vector<VertexLinkage>&& verticesLinkage);

    /**
     * Builds trees from set of possibly nested and touching, but not intersecting polygons
     * @param vertices vertex coordinates
     * @param polygons input polygons, consumed by this function, see decomposePolygonGraph()
     * @return root polygons of resulting trees
     */
    std::vector<PolygonTree> buildPolygonTrees(const std::vector<glm::dvec2>& vertices,
                                               std::vector<Polygon>&& polygons);

    /**
     * Rebuilds polygon trees to represent actual hierarchy of polygon areas rather than hierarchy of its contours
     * @param tree input polygon tree, consumed by this function, see buildPolygonTrees()
     * @return root polygons of resulting trees
     */
    std::vector<PolygonWithHolesTree> buildPolygonAreaTrees(std::vector<PolygonTree>&& tree);

    /**
     * Optimally triangulates given polygon. Optimally means that it will try to minimize number of triangles and
     * keep all triangles as close to equilateral as possible
     * @param vertices vertex coordinates
     * @param polygon input polygon, see buildPolygonAreaTrees()
     * @return vertex indices of generated triangles. Winding order of triangles is derived from parent polygon.
     * Number of triangles will not be more than 'vertices.size() + polygon.holes.size() * 2 - 2'
     */
    std::vector<glm::ivec3> triangulatePolygonWithHoles(const std::vector<glm::dvec2>& vertices,
                                                    const PolygonWithHoles& polygon);


}