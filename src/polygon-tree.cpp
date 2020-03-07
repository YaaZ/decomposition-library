#include <stdlib.h>

#include "decomposition-base.h"


namespace decomposition {


    namespace internal {


        static const double DOUBLE_EPSILON = 1E-6;


        /// This function returns 1 if edge is inside polygon, 0 if it's on boundary and -1 if it's outside
        inline int calculateEdgeToPolygonNestingByDirections(glm::dvec2 previousPolygonVertex,
                                                             glm::dvec2 currentPolygonVertex,
                                                             glm::dvec2 nextPolygonVertex,
                                                             WindingOrder polygonWinding,
                                                             glm::dvec2 anotherEdgePoint) {
            glm::dvec2 polygonEdge1 = currentPolygonVertex - previousPolygonVertex, polygonEdge2 = nextPolygonVertex - currentPolygonVertex;
            if(glm::abs(polygonEdge1.x) < DOUBLE_EPSILON && glm::abs(polygonEdge1.y) < DOUBLE_EPSILON) return 0;
            if(glm::abs(polygonEdge2.x) < DOUBLE_EPSILON && glm::abs(polygonEdge2.y) < DOUBLE_EPSILON) return 0;
            glm::dvec2 polygonInnerNormal1 = glm::dvec2(polygonEdge1.y, -polygonEdge1.x) * (double) polygonWinding;
            glm::dvec2 polygonInnerNormal2 = glm::dvec2(polygonEdge2.y, -polygonEdge2.x) * (double) polygonWinding;
            double distance1 = dot(anotherEdgePoint - currentPolygonVertex, polygonInnerNormal1);
            double distance2 = dot(anotherEdgePoint - currentPolygonVertex, polygonInnerNormal2);
            if(dot(polygonInnerNormal1, polygonEdge2) >= 0) { // Polygon vertex is convex
                if(distance1 > DOUBLE_EPSILON && distance2 > DOUBLE_EPSILON) return 1;
                else if(distance1 < -DOUBLE_EPSILON || distance2 < -DOUBLE_EPSILON) return -1;
                else return 0;
            }
            else { // Polygon vertex is concave
                if(distance1 > DOUBLE_EPSILON || distance2 > DOUBLE_EPSILON) return 1;
                else if(distance1 < -DOUBLE_EPSILON && distance2 < -DOUBLE_EPSILON) return -1;
                else return 0;
            }
        }
        /// This function returns 1 if point is inside polygon, 0 if it's on boundary and -1 if it's outside
        inline int calculatePointToPolygonNesting(const std::vector<glm::dvec2>& vertices,
                                                  const std::vector<int>& polygonVertexIndices, glm::dvec2 point) {
            int winding = 0;
            glm::dvec2 lastVertex = vertices[polygonVertexIndices.back()];
            for (int i = 0; i < polygonVertexIndices.size(); i++) {
                glm::dvec2 vertex = vertices[polygonVertexIndices[i]];
                bool crossingUpwards = vertex.y <= point.y;
                if(crossingUpwards == (lastVertex.y > point.y)) {
                    glm::dvec2 edge = vertex - lastVertex;
                    glm::dvec2 clockwiseNormal = glm::dvec2(edge.y, -edge.x);
                    double signedDistance = dot(point - vertex, clockwiseNormal);
                    if(glm::abs(signedDistance) <= DOUBLE_EPSILON) return 0;
                    else if(crossingUpwards == (signedDistance > 0)) winding += crossingUpwards ? 1 : -1;
                }
                else if(point == vertex || (vertex.y == point.y && lastVertex.y == point.y &&
                        glm::min(lastVertex.x, vertex.x) <= point.x && glm::max(lastVertex.x, vertex.x) >= point.x)) return 0;
                lastVertex = vertex;
            }
            return winding == 0 ? -1 : 1;
        }

        static bool arePolygonsNested(const std::vector<glm::dvec2>& vertices, const Polygon& parent, const Polygon& child) {
            // Most simple case: check AABB nesting
            if(
                    child.aabb.min.x < parent.aabb.min.x - DOUBLE_EPSILON || child.aabb.min.y < parent.aabb.min.y - DOUBLE_EPSILON ||
                    child.aabb.max.x > parent.aabb.max.x + DOUBLE_EPSILON || child.aabb.max.y > parent.aabb.max.y + DOUBLE_EPSILON
                    ) return false;

            // Most common case: check if at least one of child's vertices is inside parent polygon
            for (int i = 0; i < child.vertexIndices.size(); i++) {
                int vertexIndex = child.vertexIndices[i];
                glm::dvec2 vertex = vertices[vertexIndex];
                int vertexNesting = calculatePointToPolygonNesting(vertices, parent.vertexIndices, vertex);
                if(vertexNesting != 0) return vertexNesting == 1;
            }

            // Rare case: all child's vertices belongs to parent polygon as well, check edge directions
            if(parent.vertexIndices.size() < 3) return false;
            glm::dvec2 lastChildVertex = vertices[child.vertexIndices.back()];
            for (int i = 0; i < child.vertexIndices.size(); i++) {
                int currentChildVertexIndex = child.vertexIndices[i];
                glm::dvec2 currentChildVertex = vertices[currentChildVertexIndex];
                for (int j = 0; j < parent.vertexIndices.size(); j++) {
                    if(currentChildVertexIndex == parent.vertexIndices[j]) {
                        glm::dvec2 previousParentVertex = vertices[parent.vertexIndices[(j + parent.vertexIndices.size() - 1) % parent.vertexIndices.size()]];
                        glm::dvec2 currentParentVertex = vertices[parent.vertexIndices[j]];
                        glm::dvec2 nextParentVertex = vertices[parent.vertexIndices[(j + 1) % parent.vertexIndices.size()]];
                        int edgeNesting = calculateEdgeToPolygonNestingByDirections(previousParentVertex,
                                                                                    currentParentVertex, nextParentVertex, parent.getWinding(), lastChildVertex);
                        if(edgeNesting != 0) return edgeNesting == 1;
                    }
                }
                lastChildVertex = currentChildVertex;
            }

            // Very rare case: we cannot determine polygons nesting
            return true; // Why not?
        }

        static PolygonTree* findSubtreeForPolygon(const std::vector<glm::dvec2>& vertices, PolygonTree& subtree, const Polygon& polygon) {
            if(!arePolygonsNested(vertices, subtree, polygon)) return nullptr;
            for(PolygonTree& child : subtree.childrenPolygons) {
                PolygonTree* childSubtree = findSubtreeForPolygon(vertices, child, polygon);
                if(childSubtree != nullptr) return childSubtree;
            }
            return &subtree;
        }

        static void appendPolygonToTree(const std::vector<glm::dvec2>& vertices, std::vector<PolygonTree>& trees, Polygon&& polygon) {
            WindingOrder winding = polygon.getWinding();
            PolygonTree newPolygon {std::move(polygon)};
            for(PolygonTree& tree : trees) {
                PolygonTree* subtreeToAdd = findSubtreeForPolygon(vertices, tree, newPolygon);
                if(subtreeToAdd != nullptr) {
                    newPolygon.netWinding = subtreeToAdd->netWinding + winding;
                    subtreeToAdd->childrenPolygons.push_back(std::move(newPolygon));
                    return;
                }
            }
            newPolygon.netWinding = winding;
            trees.push_back(std::move(newPolygon));
        }


    }


    DECOMPOSITION_API std::vector<PolygonTree> buildPolygonTrees(const std::vector<glm::dvec2>& vertices,
                                                                std::vector<Polygon>&& polygons) {
        std::vector<PolygonTree> trees;
        std::qsort(polygons.data(), polygons.size(), sizeof(Polygon), [](const void* a, const void* b) {
            double v1 = ((Polygon*) a)->getArea(), v2 = ((Polygon*) b)->getArea();
            return (v1 < v2) ? -1 : (v1 > v2);
        });
        for (int i = polygons.size() - 1; i >= 0 ; i--) {
            internal::appendPolygonToTree(vertices, trees, std::move(polygons[i]));
        }
        return trees;
    }





    namespace internal {


        static PolygonWithHolesTree transformPolygonTreeIntoPolygonWithHolesTree(
                std::vector<PolygonTree>& siblings, PolygonTree&& tree) {
            std::vector<PolygonTree> childrenPolygons = std::move(tree.childrenPolygons);
            PolygonWithHolesTree resultTree {std::move(tree)};
            resultTree.netWinding = tree.netWinding;
            while(!childrenPolygons.empty()) {
                PolygonTree child = std::move(childrenPolygons.back());
                childrenPolygons.pop_back();
                if(child.getWinding() * resultTree.getWinding() >= 0) {
                    resultTree.childrenPolygons.push_back(
                            transformPolygonTreeIntoPolygonWithHolesTree(childrenPolygons, std::move(child))
                    );
                }
                else {
                    siblings.reserve(siblings.size() + child.childrenPolygons.size());
                    for (int i = 0; i < child.childrenPolygons.size(); i++) {
                        siblings.push_back(std::move(child.childrenPolygons[i]));
                    }
                    child.childrenPolygons.clear();
                    resultTree.holes.push_back(std::move(child));
                }
            }
            return resultTree;
        }


    }


    DECOMPOSITION_API std::vector<PolygonWithHolesTree> buildPolygonAreaTrees(std::vector<PolygonTree>&& tree) {
        std::vector<PolygonWithHolesTree> resultTree;
        resultTree.reserve(tree.size());
        while(!tree.empty()) {
            PolygonTree child = std::move(tree.back());
            tree.pop_back();
            resultTree.push_back(
                    internal::transformPolygonTreeIntoPolygonWithHolesTree(tree, std::move(child))
            );
        }
        return resultTree;
    }


}