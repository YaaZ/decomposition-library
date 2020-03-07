#include <optional>
#include <stdlib.h>

#include "decomposition-base.h"


namespace decomposition {


    namespace internal {


        struct VisitedVertex {
            int index;
            double area;
        };

        inline double calculateTrapezoidAreaByEdge(glm::dvec2 vertex1, glm::dvec2 vertex2) {
            return (vertex2.x - vertex1.x) * (vertex1.y + vertex2.y) / 2;
        }

        static std::optional<Polygon> cutLoop(const std::vector<glm::dvec2>& vertices,
                std::vector<VisitedVertex>& visitedVertices, int lastVertexIndex) {
            int firstLoopVertex = visitedVertices.size() - 1;
            Polygon::AABB aabb {vertices[lastVertexIndex], vertices[lastVertexIndex]};
            while(firstLoopVertex >= 0 && visitedVertices[firstLoopVertex].index != lastVertexIndex) {
                aabb.expand(vertices[visitedVertices[firstLoopVertex].index]);
                firstLoopVertex--;
            }
            if(firstLoopVertex == -1) return std::nullopt; // There's no loop
            double area = calculateTrapezoidAreaByEdge(
                    vertices[visitedVertices[visitedVertices.size() - 1].index],
                    vertices[lastVertexIndex]
            ) + visitedVertices[visitedVertices.size() - 1].area - visitedVertices[firstLoopVertex].area;
            Polygon resultPolygon {std::vector<int>(visitedVertices.size() - firstLoopVertex), aabb, area};
            for (int i = 0; i < resultPolygon.vertexIndices.size(); i++) {
                resultPolygon.vertexIndices[i] = visitedVertices[firstLoopVertex + i].index;
            }
            DEBUG_LOG("Cut loop from graph: " << (visitedVertices.size() - firstLoopVertex) << " last vertices")
            visitedVertices.resize(firstLoopVertex);
            return resultPolygon;
        }


        /// Returns true if a is more clockwise than b relative to negative y
        inline bool comparePathDirections(const glm::dvec2 a, const glm::dvec2 b) {
            if(b.x >= 0 && a.x >= 0) return b.y < a.y;
            else if(b.x < 0 && a.x < 0) return b.y > a.y;
            else return b.x >= 0;
        }
        /// This function selects most left/right not yet visited path depending on passed clockwise/counter-clockwise order respectively
        inline VertexLinkage::Path& choosePath(VertexLinkage& vertex, const glm::dvec2 inDirection, const WindingOrder order) {
            assert(vertex.pathsNumber != 0);
            if(vertex.pathsNumber == 1) return vertex.getPath(0);
            /* This matrix rotates (and probably reflects) vectors so inDirection points upwards and positive x
             * is preferable direction according to order parameter */
            glm::dmat2 transformMatrix = glm::dmat2(inDirection.y * (double) order, inDirection.x, -inDirection.x * (double) order, inDirection.y);
            int bestPathIndex = -1;
            glm::dvec2 bestPathDirection(-1, -2);
            for (int i = 0; i < vertex.pathsNumber; i++) {
                VertexLinkage::Path& path = vertex.getPath(i);
                if(path.next == -1) continue;
                glm::dvec2 pathDirection = transformMatrix * path.directionToNextVertex;
                if(comparePathDirections(bestPathDirection, pathDirection)) {
                    bestPathIndex = i;
                    bestPathDirection = pathDirection;
                }
            }
            assert(bestPathIndex != -1);
            return vertex.getPath(bestPathIndex);
        }

        static void traverseGraph(const std::vector<glm::dvec2>& vertices, std::vector<Polygon>& polygons,
                std::vector<VertexLinkage>& verticesLinkage, std::vector<int>& lastVisitIterationsByVertex,
                int iteration, int startingVertexIndex, int startingPathIndex) {
            DEBUG_LOG(std::endl << "Starting graph traversing from vertex: " << startingVertexIndex)
            std::vector<VisitedVertex> visitedVertices {{startingVertexIndex, 0}};
            VertexLinkage& startingVertex = verticesLinkage[startingVertexIndex];
            lastVisitIterationsByVertex[startingVertexIndex] = iteration;
            VertexLinkage::Path& startingPath = startingVertex.getPath(startingPathIndex);
            int nextVertexIndex = startingPath.next;
            startingPath.next = -1;
            glm::dvec2 currentDirection = startingPath.directionToNextVertex;
            for(;;) {
                int vertexIndex = nextVertexIndex;
                VertexLinkage& vertex = verticesLinkage[vertexIndex];
                if(lastVisitIterationsByVertex[vertexIndex] == iteration) {
                    std::optional<Polygon> polygon = cutLoop(vertices, visitedVertices, vertexIndex);
                    if(polygon) polygons.push_back(std::move(*polygon));
                    if(visitedVertices.empty()) break;
                }
                DEBUG_ONLY(
                        std::cout << "Visited vertex: " << vertexIndex << ", paths:";
                        for (int i = 0; i < vertex.pathsNumber; i++) {
                            std::cout << " " << vertex.getPath(i).next;
                        }
                        std::cout << std::endl;
                )
                glm::dvec2 vertexCoordinates = vertices[vertexIndex];
                lastVisitIterationsByVertex[vertexIndex] = iteration;
                auto [previousVertexIndex, previousVertexArea] = visitedVertices[visitedVertices.size() - 1];
                glm::dvec2 previousVertexCoordinates = vertices[previousVertexIndex];
                visitedVertices.push_back({
                    vertexIndex,
                    previousVertexArea + calculateTrapezoidAreaByEdge(previousVertexCoordinates, vertexCoordinates)
                });
                /* We don't really care about order: it's just a hint, so it just needs to be constant while
                 * algorithm is running. So we just pass clockwise */
                VertexLinkage::Path& chosenPath = choosePath(vertex, currentDirection, WindingOrder::CLOCKWISE);
                nextVertexIndex = chosenPath.next;
                chosenPath.next = -1;
                currentDirection = chosenPath.directionToNextVertex;
            }
        }


        struct SortableVertex {
            double x;
            int index;
        };

        static std::vector<SortableVertex> sortVertices(const std::vector<glm::dvec2>& vertices) {
            std::vector<SortableVertex> sortedVertices(vertices.size());
            for (int i = 0; i < vertices.size(); i++) {
                sortedVertices[i] = {vertices[i].x, i};
            }
            std::qsort(sortedVertices.data(), sortedVertices.size(), sizeof(SortableVertex), [](const void* a, const void* b) {
                double v1 = ((SortableVertex*) a)->x, v2 = ((SortableVertex*) b)->x;
                return (v1 < v2) ? -1 : (v1 > v2);
            });
            return sortedVertices;
        }


    }


    DECOMPOSITION_API std::vector<Polygon> decomposePolygonGraph(const std::vector<glm::dvec2>& vertices,
                                                                 std::vector<VertexLinkage>&& verticesLinkage) {
        DEBUG_LOG(std::endl << "Starting graph decomposition")
        std::vector<Polygon> polygons;
        std::vector<int> lastVisitIterationsByVertex(vertices.size(), -1);
        int iteration = 0;
        std::vector<internal::SortableVertex> sortedVertices = internal::sortVertices(vertices);
        for (int i = 0; i < sortedVertices.size(); i++) {
            int vertexIndex = sortedVertices[i].index;
            VertexLinkage& linkage = verticesLinkage[vertexIndex];
            for (int j = 0; j < linkage.pathsNumber; j++) {
                if(linkage.getPath(j).next != -1) {
                    internal::traverseGraph(vertices, polygons, verticesLinkage, lastVisitIterationsByVertex, iteration++, vertexIndex, j);
                }
            }
        }
        return polygons;
    }


}