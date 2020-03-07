#include "decomposition-base.h"


namespace decomposition {


    namespace internal {


        static const double DOUBLE_EPSILON = 1E-6;


        struct Vertex {

            int index, previous, next;
            glm::dvec2 relativeDirection;
            std::vector<int> cuts;

            inline bool isConvex() const {
                return relativeDirection.x > 0;
            }

        };

        struct VertexGraph {

            struct Polygon {
                int mostRightVertexIndex, indexOffsetInVertices;
                decomposition::Polygon::AABB aabb;
            };

            std::vector<Vertex> vertices;
            std::vector<Polygon> polygonsInfo;

        };


        /** This matrix rotates (and probably reflects) vectors so inDirection points upwards and positive x
         *  is convex angle and negative x is concave */
        inline glm::dmat2 createDirectionTransformationMatrix(const glm::dvec2& inDirection, const WindingOrder winding) {
            /* This matrix rotates (and probably reflects) vectors so inDirection points upwards and positive x
             * is convex angle and negative x is concave */
            glm::dvec2 dir = glm::normalize(inDirection);
            return glm::dmat2(dir.y * (double) winding, dir.x, -dir.x * (double) winding, dir.y);
        }
        /// Returns true if b is more clockwise than a relative to negative y
        inline bool compareDirections(const glm::dvec2& a, const glm::dvec2& b, const bool oppositeFirst) {
            // There is a compilation bug with inlined function when passing arguments by value
            if(glm::abs(a.x) <= DOUBLE_EPSILON && a.y < 0) return !oppositeFirst;
            else if(glm::abs(b.x) <= DOUBLE_EPSILON && b.y < 0) return oppositeFirst;
            else if(b.x > 0 && a.x > 0) return b.y < a.y;
            else if(b.x < 0 && a.x < 0) return b.y > a.y;
            else return b.x > a.x;
        }





        static int buildVertexGraph(const std::vector<glm::dvec2>& vertices, std::vector<Vertex>& outputVertices,
                                    const Polygon& polygon, const WindingOrder winding) {
            // Put vertices into output list, deleting repeating vertices
            glm::dvec2 previousVertex = vertices[polygon.vertexIndices.back()];
            int firstOutputVertexIndex = outputVertices.size();
            int mostRightVertexIndex = -1;
            double mostRightVertexCoordinate = previousVertex.x;
            for (int i = 0; i < polygon.vertexIndices.size(); i++) {
                int vertexIndex = polygon.vertexIndices[i];
                glm::dvec2 vertex = vertices[vertexIndex];
                glm::dvec2 edge = vertex - previousVertex;
                if(i == 0 || glm::abs(edge.x) > DOUBLE_EPSILON || glm::abs(edge.y) > DOUBLE_EPSILON) {
                    if(i == 0 || vertex.x > mostRightVertexCoordinate) {
                        mostRightVertexCoordinate = vertex.x;
                        mostRightVertexIndex = outputVertices.size();
                    }
                    outputVertices.push_back({vertexIndex, -1, (int) outputVertices.size() + 1});
                }
                previousVertex = vertex;
            }
            outputVertices.back().next = firstOutputVertexIndex;

            // Calculate convexity
            int addedVertices = outputVertices.size() - firstOutputVertexIndex;
            previousVertex = vertices[outputVertices.back().index];
            glm::dvec2 currentVertex = vertices[outputVertices[firstOutputVertexIndex].index];
            for (int i = 0; i < addedVertices; i++) {
                glm::dvec2 nextVertex = vertices[outputVertices[firstOutputVertexIndex + (i + 1) % addedVertices].index];
                outputVertices[firstOutputVertexIndex + i].relativeDirection =
                        createDirectionTransformationMatrix(currentVertex - previousVertex, winding) *
                        glm::normalize(nextVertex - currentVertex);
                previousVertex = currentVertex;
                currentVertex = nextVertex;
            }
            return mostRightVertexIndex;
        }

        static VertexGraph buildVertexGraph(const std::vector<glm::dvec2>& vertices, const PolygonWithHoles& polygon) {
            int totalVertexNumber = polygon.vertexIndices.size();
            for(const Polygon& hole : polygon.holes) totalVertexNumber += hole.vertexIndices.size();
            VertexGraph graph {};
            graph.vertices.reserve(totalVertexNumber + polygon.holes.size() * 2);
            graph.polygonsInfo.reserve(polygon.holes.size() + 1);
            for(const Polygon& hole : polygon.holes) {
                int offset = graph.vertices.size();
                int mostRightVertexIndex = buildVertexGraph(vertices, graph.vertices, hole, polygon.getWinding());
                if(mostRightVertexIndex != -1) graph.polygonsInfo.push_back({mostRightVertexIndex, offset, hole.aabb});
            }
            int offset = graph.vertices.size();
            int mostRightVertexIndex = buildVertexGraph(vertices, graph.vertices, polygon, polygon.getWinding());
            if(mostRightVertexIndex != -1) graph.polygonsInfo.push_back({mostRightVertexIndex, offset, polygon.aabb});
            else graph.polygonsInfo.clear();
            return graph;
        }





        struct MutualVisibilityCandidate {
            glm::dvec2 lookupVertex, intersectionPoint, candidateVertex, rejectedCandidateVertex;
            int lookupVertexIndex, candidateVertexIndex, rejectedCandidateVertexIndex;
        };

        static MutualVisibilityCandidate findVisibilityCandidate(const std::vector<glm::dvec2>& vertices, const VertexGraph& graph,
                                                                 const int polygonIndex, const WindingOrder winding) {
            int lookupVertexIndex = graph.polygonsInfo[polygonIndex].mostRightVertexIndex;
            DEBUG_LOG(std::endl << "Lookup vertex:  " << graph.vertices[lookupVertexIndex].index)
            glm::dvec2 lookupVertex = vertices[graph.vertices[lookupVertexIndex].index];
            glm::dvec2 closestRightIntersectionPoint(HUGE_VAL, lookupVertex.y);
            int closestRightIntersectionPointEdge = -1;
            for (int i = 0; i < graph.polygonsInfo.size(); i++) {
                if(i == polygonIndex) continue;
                Polygon::AABB aabb = graph.polygonsInfo[i].aabb;
                if(aabb.max.x < lookupVertex.x || aabb.max.y < lookupVertex.y || aabb.min.y > lookupVertex.y) continue;
                int firstVertexIndex = graph.polygonsInfo[i].indexOffsetInVertices;
                glm::dvec2 lastVertex = vertices[graph.vertices[firstVertexIndex].index];
                int vertexIndex = firstVertexIndex;
                do {
                    const Vertex& vertex = graph.vertices[vertexIndex];
                    glm::dvec2 currentVertex = lastVertex;
                    glm::dvec2 nextVertex = vertices[graph.vertices[vertex.next].index];
                    if((currentVertex.x >= lookupVertex.x || nextVertex.x >= lookupVertex.x) &&
                       (nextVertex.y - lookupVertex.y) * ((double) winding) <= 0 &&
                       (currentVertex.y - lookupVertex.y) * (nextVertex.y - lookupVertex.y) <= 0) {
                        double intersectionPointOffset = (currentVertex.y - lookupVertex.y) / (currentVertex.y - nextVertex.y);
                        glm::dvec2 intersectionPoint = currentVertex + (nextVertex - currentVertex) * intersectionPointOffset;
                        DEBUG_LOG("  Candidate edge: " << vertex.index << " " << graph.vertices[vertex.next].index)
                        if(intersectionPoint.x >= lookupVertex.x && intersectionPoint.x < closestRightIntersectionPoint.x) {
                            closestRightIntersectionPoint = intersectionPoint;
                            closestRightIntersectionPointEdge = vertexIndex;
                        }
                    }
                    lastVertex = nextVertex;
                    vertexIndex = vertex.next;
                } while(vertexIndex != firstVertexIndex);
            }
            assert(closestRightIntersectionPointEdge != -1);
            int candidateVertexIndex1 = closestRightIntersectionPointEdge;
            int candidateVertexIndex2 = graph.vertices[closestRightIntersectionPointEdge].next;
            glm::dvec2 candidateVertex1 = vertices[graph.vertices[candidateVertexIndex1].index];
            glm::dvec2 candidateVertex2 = vertices[graph.vertices[candidateVertexIndex2].index];
            DEBUG_LOG("  Closest edge: " << graph.vertices[closestRightIntersectionPointEdge].index << " "
                      << graph.vertices[graph.vertices[closestRightIntersectionPointEdge].next].index)
            bool pickCandidate1 = candidateVertex1.x >= candidateVertex2.x;
            DEBUG_LOG("  Visibility candidate: " << graph.vertices[pickCandidate1 ? candidateVertexIndex1 : candidateVertexIndex2].index)
            return {
                    lookupVertex,
                    closestRightIntersectionPoint,
                    pickCandidate1 ? candidateVertex1 : candidateVertex2,
                    pickCandidate1 ? candidateVertex2 : candidateVertex1,
                    lookupVertexIndex,
                    pickCandidate1 ? candidateVertexIndex1 : candidateVertexIndex2,
                    pickCandidate1 ? candidateVertexIndex2 : candidateVertexIndex1
            };
        }

        static int findVisibleVertex(const std::vector<glm::dvec2>& vertices, const VertexGraph& graph,
                                     const int polygonIndex, const MutualVisibilityCandidate& candidate) {
            if(candidate.candidateVertex.y == candidate.lookupVertex.y) return candidate.candidateVertexIndex;
            double verticalInnerDirection = candidate.candidateVertex.y > candidate.lookupVertex.y ? 1 : -1;
            glm::dvec2 vectorToCandidate = candidate.candidateVertex - candidate.lookupVertex;
            glm::dvec2 lookupInnerNormal = glm::dvec2(vectorToCandidate.y, -vectorToCandidate.x) * verticalInnerDirection;
            glm::dvec2 candidateEdgeDirection = candidate.rejectedCandidateVertex - candidate.candidateVertex;
            glm::dvec2 candidateInnerNormal = glm::dvec2(candidateEdgeDirection.y, -candidateEdgeDirection.x) * verticalInnerDirection;
            Polygon::AABB triangleAabb = {
                    glm::min(glm::min(candidate.lookupVertex, candidate.intersectionPoint), candidate.candidateVertex),
                    glm::max(glm::max(candidate.lookupVertex, candidate.intersectionPoint), candidate.candidateVertex)
            };
            int closestVertexIndex = -1;
            double closestVertexTg = 0, closestVertexDistance = 0;
            for (int i = 0; i < graph.polygonsInfo.size(); i++) {
                if(i == polygonIndex) continue;
                Polygon::AABB aabb = graph.polygonsInfo[i].aabb;
                if(aabb.max.x < triangleAabb.min.x - DOUBLE_EPSILON || aabb.max.y < triangleAabb.min.y - DOUBLE_EPSILON ||
                   aabb.min.x > triangleAabb.max.x + DOUBLE_EPSILON || aabb.min.y > triangleAabb.max.y + DOUBLE_EPSILON) continue;
                int firstVertexIndex = graph.polygonsInfo[i].indexOffsetInVertices;
                int lastVertexIndex = i == graph.polygonsInfo.size() - 1 ? graph.vertices.size() : graph.polygonsInfo[i + 1].indexOffsetInVertices;
                for (int j = firstVertexIndex; j < lastVertexIndex; j++) {
                    const Vertex& vertex = graph.vertices[j];
                    if(vertex.isConvex()) continue;
                    glm::dvec2 v = vertices[vertex.index];
                    glm::dvec2 vertexFromLookupPoint = v - candidate.lookupVertex;
                    double verticalOffset = vertexFromLookupPoint.y * verticalInnerDirection;
                    if(verticalOffset >= -DOUBLE_EPSILON &&
                       glm::dot(vertexFromLookupPoint, lookupInnerNormal) >= -DOUBLE_EPSILON &&
                       glm::dot(v - candidate.candidateVertex, candidateInnerNormal) >= -DOUBLE_EPSILON) {
                        DEBUG_LOG("  New candidate: " << vertex.index)
                        double vertexTg = verticalOffset / (v.x - candidate.lookupVertex.x);
                        double vertexDistance = glm::length(v - candidate.lookupVertex);
                        if(closestVertexIndex == -1 || vertexTg < closestVertexTg - DOUBLE_EPSILON ||
                           (glm::abs(vertexTg - closestVertexTg) <= DOUBLE_EPSILON &&
                           vertexDistance < closestVertexDistance)) {
                            closestVertexIndex = j;
                            closestVertexTg = vertexTg;
                            closestVertexDistance = vertexDistance;
                        }
                    }
                }
            }
            return closestVertexIndex == -1 ? candidate.candidateVertexIndex : closestVertexIndex;
        }

        inline void cutGraph(VertexGraph& graph, int vertex1, int vertex2) {
            if(graph.vertices[vertex1].next != -1) {
                graph.vertices[vertex1].cuts.push_back(graph.vertices[vertex1].next);
                graph.vertices[vertex1].next = -1;
            }
            if(graph.vertices[vertex2].next != -1) {
                graph.vertices[vertex2].cuts.push_back(graph.vertices[vertex2].next);
                graph.vertices[vertex2].next = -1;
            }
            graph.vertices[vertex1].cuts.push_back(vertex2);
            graph.vertices[vertex2].cuts.push_back(vertex1);
        }
        static void cutPolygonGraphIntoSimple(const std::vector<glm::dvec2>& vertices, VertexGraph& graph, const WindingOrder winding) {
            std::vector<glm::ivec2> cuts(graph.polygonsInfo.size() - 1);
            for (int i = 0; i < graph.polygonsInfo.size() - 1; i++) {
                MutualVisibilityCandidate candidate = findVisibilityCandidate(vertices, graph, i, winding);
                int visibleVertexIndex;
                if(glm::abs(candidate.lookupVertex.y - candidate.candidateVertex.y) <= DOUBLE_EPSILON) {
                    visibleVertexIndex = candidate.candidateVertexIndex;
                }
                else if(glm::abs(candidate.lookupVertex.y - candidate.rejectedCandidateVertex.y) <= DOUBLE_EPSILON) {
                    visibleVertexIndex = candidate.rejectedCandidateVertexIndex;
                }
                else visibleVertexIndex = findVisibleVertex(vertices, graph, i, candidate);
                DEBUG_LOG("Visible vertex: " << graph.vertices[visibleVertexIndex].index)
                cuts[i] = {candidate.lookupVertexIndex, visibleVertexIndex};
            }
            for (int i = 0; i < cuts.size(); i++) cutGraph(graph, cuts[i].x, cuts[i].y);
        }





        /// This function selects most right/left not yet visited path depending on passed clockwise/counter-clockwise order respectively
        inline int choosePath(const std::vector<glm::dvec2>& vertices, int previousVertexIndex, int vertexIndex, VertexGraph& graph, const WindingOrder winding) {
            Vertex& vertex = graph.vertices[vertexIndex];
            glm::dvec2 currentVertex = vertices[vertex.index];
            glm::dmat2 directionMatrix = createDirectionTransformationMatrix(
                    currentVertex - vertices[graph.vertices[previousVertexIndex].index], winding);
            int bestPathIndex = -1;
            glm::dvec2 bestPathDirection;
            for (int i = 0; i < vertex.cuts.size(); i++) {
                Vertex& cut = graph.vertices[vertex.cuts[i]];
                DEBUG_LOG("  Possible path: " << cut.index)
                glm::dvec2 cutVertex = vertices[cut.index];
                glm::dvec2 path = cutVertex - currentVertex;
                if(glm::abs(path.x) <= DOUBLE_EPSILON && glm::abs(path.y) <= DOUBLE_EPSILON) {
                    DEBUG_LOG("    Merged cut vertex")
                    if(cut.cuts.size() != 0) {
                        for (int c : cut.cuts) {
                            if(c != vertexIndex) vertex.cuts.push_back(c);
                        }
                        cut.cuts.clear();
                    }
                    else if(cut.next != vertexIndex) vertex.cuts.push_back(cut.next);
                    cut.next = vertexIndex;
                    vertex.cuts.erase(vertex.cuts.begin() + i);
                    i--;
                    continue;
                }
                glm::dvec2 pathDirection = directionMatrix * glm::normalize(path);
                if(bestPathIndex == -1 || compareDirections(bestPathDirection, pathDirection, false)) {
                    bestPathIndex = i;
                    bestPathDirection = pathDirection;
                }
            }
            assert(bestPathIndex != -1);
            return bestPathIndex;
        }

        static void buildSimplePolygon(const std::vector<glm::dvec2>& vertices,
                                       VertexGraph& graph, const WindingOrder winding) {
            int startingVertexIndex = graph.polygonsInfo.back().indexOffsetInVertices;
            while(!graph.vertices[startingVertexIndex].cuts.empty()) startingVertexIndex++;
            DEBUG_LOG(std::endl << "Starting graph traversing from vertex: " <<
                      graph.vertices[startingVertexIndex].index << ", " <<
                      (graph.vertices[startingVertexIndex].isConvex() ? "convex" : "concave"))
            int currentVertexIndex = startingVertexIndex;
            bool traversedFullCycle = false;
            for(;;) {
                int previousVertexIndex = currentVertexIndex;
                Vertex& previousVertex = graph.vertices[previousVertexIndex];
                currentVertexIndex = previousVertex.next;
                Vertex& currentVertex = graph.vertices[currentVertexIndex];
                if(currentVertexIndex == startingVertexIndex) traversedFullCycle = true;
                DEBUG_LOG("Visited vertex: " << currentVertex.index << ", " <<
                    (currentVertex.cuts.size() == 0 ? (currentVertex.isConvex() ? "convex" : "concave") : "cutting point"))
                if(currentVertex.cuts.size() == 0) {
                    glm::dvec2 nextEdge = vertices[graph.vertices[graph.vertices[currentVertexIndex].next].index] -
                                          vertices[graph.vertices[currentVertexIndex].index];
                    if(glm::abs(nextEdge.x) < DOUBLE_EPSILON && glm::abs(nextEdge.y) < DOUBLE_EPSILON) {
                        if(graph.vertices[currentVertexIndex].next == startingVertexIndex) traversedFullCycle = true;
                        Vertex& newCurrentVertex = graph.vertices[graph.vertices[currentVertexIndex].next];
                        previousVertex.next = graph.vertices[currentVertexIndex].next;
                        if(newCurrentVertex.next != -1) {
                            glm::dmat2 directionMatrix = createDirectionTransformationMatrix(
                                    vertices[newCurrentVertex.index] - vertices[previousVertex.index], winding);
                            newCurrentVertex.relativeDirection = directionMatrix *
                                                                 glm::normalize(vertices[graph.vertices[newCurrentVertex.next].index] - vertices[newCurrentVertex.index]);
                        }
                        DEBUG_LOG("  Deleted repeating vertex")
                        if(previousVertex.next == newCurrentVertex.next) break;
                        currentVertexIndex = previousVertexIndex;
                    }
                    else if(traversedFullCycle) break;
                }
                else {
                    int cutIndex = choosePath(vertices, previousVertexIndex, currentVertexIndex, graph, winding);
                    int cut = currentVertex.cuts[cutIndex];
                    currentVertex.cuts.erase(currentVertex.cuts.begin() + cutIndex, currentVertex.cuts.begin() + cutIndex + 1);

                    glm::dmat2 directionMatrix = createDirectionTransformationMatrix(
                            vertices[currentVertex.index] - vertices[previousVertex.index], winding);
                    glm::dvec2 relativeDirection = directionMatrix *
                            glm::normalize(vertices[graph.vertices[cut].index] - vertices[currentVertex.index]);

                    if(currentVertex.cuts.size() == 0) {
                        currentVertex.relativeDirection = relativeDirection;
                        currentVertex.next = cut;
                    }
                    else {
                        int newVertexIndex = graph.vertices.size();
                        graph.vertices.push_back({currentVertex.index, -1, cut});
                        graph.vertices.back().relativeDirection = relativeDirection;
                        currentVertexIndex = previousVertex.next = newVertexIndex;
                    }
                    DEBUG_LOG("  New generated edge: " << currentVertex.index << " - " << graph.vertices[cut].index <<
                                                       ", " << (graph.vertices.back().isConvex() ? "convex" : "concave"))
                }
            }
        }





        static bool isEar(const std::vector<glm::dvec2>& vertices, const VertexGraph& graph, const WindingOrder winding,
                          const Vertex& currentVertex, bool convex) {
            const Vertex& previousVertex = graph.vertices[currentVertex.previous];
            const Vertex& nextVertex = graph.vertices[currentVertex.next];
            glm::dvec2 vertex1 = vertices[previousVertex.index];
            glm::dvec2 vertex2 = vertices[currentVertex.index];
            glm::dvec2 vertex3 = vertices[nextVertex.index];
            glm::dvec2 edge12 = vertex2 - vertex1, edge23 = vertex3 - vertex2, edge31 = vertex1 - vertex3;
            glm::dvec2 innerNormal12 = glm::dvec2(edge12.y, -edge12.x) * (double) winding;
            glm::dvec2 innerNormal23 = glm::dvec2(edge23.y, -edge23.x) * (double) winding;
            glm::dvec2 innerNormal31 = glm::dvec2(edge31.y, -edge31.x) * (double) winding;
            int vertexIndex = nextVertex.next;
            for(;;) {
                if(vertexIndex == currentVertex.previous) return true;
                const Vertex& vertex = graph.vertices[vertexIndex];
                if(vertex.isConvex() != convex) {
                    glm::dvec2 v = vertices[vertex.index];
                    if(glm::dot(v - vertex1, innerNormal12) >= DOUBLE_EPSILON &&
                       glm::dot(v - vertex2, innerNormal23) >= DOUBLE_EPSILON &&
                       glm::dot(v - vertex3, innerNormal31) >= DOUBLE_EPSILON) return false;
                }
                vertexIndex = vertex.next;
            }
        }

        static void cutTriangle(const std::vector<glm::dvec2>& vertices, VertexGraph& graph,
                                std::vector<glm::ivec3>& result, Vertex& currentVertex, const WindingOrder winding) {
            Vertex& previousVertex = graph.vertices[currentVertex.previous];
            Vertex& nextVertex = graph.vertices[currentVertex.next];
            DEBUG_LOG("Cut triangle: " << previousVertex.index << " " << currentVertex.index << " " << nextVertex.index)
            result.emplace_back(previousVertex.index, currentVertex.index, nextVertex.index);
            previousVertex.next = currentVertex.next;
            nextVertex.previous = currentVertex.previous;
            glm::dvec2 previousEdge = vertices[previousVertex.index] - vertices[graph.vertices[previousVertex.previous].index];
            glm::dvec2 currentEdge = vertices[nextVertex.index] - vertices[previousVertex.index];
            glm::dvec2 nextEdge = vertices[graph.vertices[nextVertex.next].index] - vertices[nextVertex.index];
            previousVertex.relativeDirection =
                    createDirectionTransformationMatrix(previousEdge, winding) * currentEdge;
            nextVertex.relativeDirection =
                    createDirectionTransformationMatrix(currentEdge, winding) * nextEdge;
        }

        inline int binarySearchVertexByDirection(const std::vector<int>& sortedVertices,
                                                 const std::vector<Vertex>& vertices,
                                                 const glm::dvec2& targetDirection) {
            int left = 0, right = sortedVertices.size() - 1;
            int center = 0;
            glm::dvec2 direction = targetDirection;
            while(left <= right) {
                center = (left + right) / 2;
                direction = vertices[sortedVertices[center]].relativeDirection;
                if(compareDirections(direction, targetDirection, true)) left = center + 1;
                else right = center - 1;
            }
            return compareDirections(direction, targetDirection, true) ? center + 1 : center;
        }
        static std::vector<int> sortVerticesByDirections(VertexGraph& graph) {
            std::vector<int> sortedVertices;
            sortedVertices.reserve(graph.vertices.size() - 2);
            int startingVertexIndex = graph.polygonsInfo.back().indexOffsetInVertices;
            int previousVertexIndex = startingVertexIndex;
            int currentVertexIndex = graph.vertices[startingVertexIndex].next;
            do {
                Vertex& currentVertex = graph.vertices[currentVertexIndex];
                currentVertex.previous = previousVertexIndex;
                sortedVertices.insert(
                        sortedVertices.begin() +
                        binarySearchVertexByDirection(sortedVertices, graph.vertices, currentVertex.relativeDirection),
                        currentVertexIndex);
                previousVertexIndex = currentVertexIndex;
                currentVertexIndex = currentVertex.next;
            } while(previousVertexIndex != startingVertexIndex);
            return sortedVertices;
        }

        static std::vector<glm::ivec3> triangulatePolygonGraph(const std::vector<glm::dvec2>& vertices,
                                                               VertexGraph& graph, const WindingOrder winding) {
            if(graph.vertices.size() < 3) return {};
            std::vector<int> sortedVertices = sortVerticesByDirections(graph);
            std::vector<glm::ivec3> result;
            result.reserve(graph.vertices.size() - 2);
            while(sortedVertices.size() >= 3) {
                int ear = -1;
                for (int i = sortedVertices.size() - 1; i >= 0 ; i--) {
                    if(isEar(vertices, graph, winding, graph.vertices[sortedVertices[i]], true)) {
                        ear = i;
                        break;
                    }
                }
                assert(ear != -1);
                int earIndex = sortedVertices[ear];
                sortedVertices.erase(sortedVertices.begin() + ear);
                Vertex& earVertex = graph.vertices[earIndex];
                int deletedVerticesFromSorted = 0;
                for (int i = sortedVertices.size() - 1; i >= 0 ; i--) {
                    if(sortedVertices[i] == earVertex.previous || sortedVertices[i] == earVertex.next) {
                        deletedVerticesFromSorted++;
                        sortedVertices.erase(sortedVertices.begin() + i);
                        if(deletedVerticesFromSorted >= 2) break;
                    }
                }
                cutTriangle(vertices, graph, result, earVertex, winding);
                sortedVertices.insert(
                        sortedVertices.begin() +
                        binarySearchVertexByDirection(sortedVertices, graph.vertices, graph.vertices[earVertex.previous].relativeDirection),
                        earVertex.previous);
                sortedVertices.insert(
                        sortedVertices.begin() +
                        binarySearchVertexByDirection(sortedVertices, graph.vertices, graph.vertices[earVertex.next].relativeDirection),
                        earVertex.next);
            }
            return result;
        }


    }


    DECOMPOSITION_API std::vector<glm::ivec3> triangulatePolygonWithHoles(const std::vector<glm::dvec2>& vertices,
                                                                          const PolygonWithHoles& polygon) {
        internal::VertexGraph graph = internal::buildVertexGraph(vertices, polygon);
        DEBUG_LOG(std::endl << "Cutting polygon into simple")
        internal::cutPolygonGraphIntoSimple(vertices, graph, polygon.getWinding());
        internal::buildSimplePolygon(vertices, graph, polygon.getWinding());
        DEBUG_LOG(std::endl << "Triangulating polygon")
        return internal::triangulatePolygonGraph(vertices, graph, polygon.getWinding());
    }


}