#include <optional>

#include "decomposition-base.h"


namespace decomposition {


    namespace internal {


        struct Intersection {
            glm::dvec2 point;
            double edge1Offset, edge2Offset;
        };

        inline std::optional<Intersection> findIntersectionPoint(const glm::dvec2 edge1Vertex1,
                                                                 const glm::dvec2 edge1Vertex2,
                                                                 const glm::dvec2 edge2Vertex1,
                                                                 const glm::dvec2 edge2Vertex2) {
            glm::dvec2 edge1 = edge1Vertex2 - edge1Vertex1;
            glm::dvec2 edge1RotatedBy90 {edge1.y, -edge1.x};
            double edge2Vertex1RelativeSignedDistanceFromEdge1 = dot(edge1RotatedBy90, edge2Vertex1 - edge1Vertex1);
            double edge2Vertex2RelativeSignedDistanceFromEdge1 = dot(edge1RotatedBy90, edge2Vertex2 - edge1Vertex1);
            if(edge2Vertex1RelativeSignedDistanceFromEdge1 * edge2Vertex2RelativeSignedDistanceFromEdge1 > 0 ||
               edge2Vertex2RelativeSignedDistanceFromEdge1 == 0) return std::nullopt;
            double edge2Offset = edge2Vertex1RelativeSignedDistanceFromEdge1 /
                                 (edge2Vertex1RelativeSignedDistanceFromEdge1 - edge2Vertex2RelativeSignedDistanceFromEdge1);
            glm::dvec2 intersectionPoint = edge2Vertex1 + (edge2Vertex2 - edge2Vertex1) * edge2Offset;
            double edge1Length = glm::length(edge1);
            double edge1Offset = dot(edge1, intersectionPoint - edge1Vertex1) / (edge1Length * edge1Length);
            if(edge1Offset < 0 || edge1Offset >= 1) return std::nullopt;
            else return Intersection{intersectionPoint, edge1Offset, edge2Offset};
        }



        struct IntersectionPoint {
            int vertexIndex;
            double edgeOffset; // offset (interpolation) value in projection on edge in range between 0 and 1
        };

        inline int binarySearchInIntersectionPoints(const std::vector<internal::IntersectionPoint>& intersectionPoints,
                                                    double targetOffset) {
            int left = 0, right = intersectionPoints.size() - 1;
            int center = 0;
            double offset = targetOffset;
            while(left <= right) {
                center = (left + right) / 2;
                offset = intersectionPoints[center].edgeOffset;
                if(offset < targetOffset) left = center + 1;
                else if(offset > targetOffset) right = center - 1;
                else return center;
            }
            return targetOffset > offset ? center + 1 : center;
        }

        static void addEdgesIntersectionToList(std::vector<internal::IntersectionPoint>& edge1IntersectionPoints,
                                               std::vector<internal::IntersectionPoint>& edge2IntersectionPoints,
                                               std::vector<glm::dvec2>& vertices,
                                               int edge1Index, int edge1EndIndex, int edge2Index, int edge2EndIndex) {
            std::optional<internal::Intersection> intersection = internal::findIntersectionPoint(
                    vertices[edge1Index], vertices[edge1EndIndex], vertices[edge2Index], vertices[edge2EndIndex]
            );
            if(intersection) {
                int edge1IntersectionListIndex = binarySearchInIntersectionPoints(edge1IntersectionPoints, intersection->edge1Offset);
                int edge2IntersectionListIndex = binarySearchInIntersectionPoints(edge2IntersectionPoints, intersection->edge2Offset);

                bool repeatedIntersectionOnEdge1 = edge1IntersectionListIndex < edge1IntersectionPoints.size() &&
                        intersection->edge1Offset == edge1IntersectionPoints[edge1IntersectionListIndex].edgeOffset;
                bool repeatedIntersectionOnEdge2 = edge2IntersectionListIndex < edge2IntersectionPoints.size() &&
                        intersection->edge2Offset == edge2IntersectionPoints[edge2IntersectionListIndex].edgeOffset;
                if(repeatedIntersectionOnEdge1 && repeatedIntersectionOnEdge2) return;

                int vertexIndex;
                if(repeatedIntersectionOnEdge1) vertexIndex = edge1IntersectionPoints[edge1IntersectionListIndex].vertexIndex;
                else if(repeatedIntersectionOnEdge2) vertexIndex = edge2IntersectionPoints[edge2IntersectionListIndex].vertexIndex;
                else if(intersection->edge1Offset == 0) vertexIndex = edge1Index;
                else if(intersection->edge2Offset == 0) vertexIndex = edge2Index;
                else {
                    vertexIndex = vertices.size();
                    vertices.push_back(intersection->point);
                }

                if(!repeatedIntersectionOnEdge1) {
                    edge1IntersectionPoints.insert(edge1IntersectionPoints.begin() + edge1IntersectionListIndex, {vertexIndex, intersection->edge1Offset});
                }
                if(!repeatedIntersectionOnEdge2) {
                    edge2IntersectionPoints.insert(edge2IntersectionPoints.begin() + edge2IntersectionListIndex, {vertexIndex, intersection->edge2Offset});
                }
            }
        }

        static void addPolygonSelfIntersectionsToList(std::vector<glm::dvec2>& vertices, const std::vector<int>& vertexIndices,
                                                      std::vector<std::vector<internal::IntersectionPoint>>& intersectionPointsByEdge) {
            for (int i = 0; i < vertexIndices.size() - 2; i++) {
                std::vector<internal::IntersectionPoint>& edge1IntersectionPoints = intersectionPointsByEdge[i];
                int checkIntersectionsTill = i == 0 ? vertexIndices.size() - 1 : vertexIndices.size();
                for (int j = i + 2; j < checkIntersectionsTill; j++) {
                    std::vector<internal::IntersectionPoint>& edge2IntersectionPoints = intersectionPointsByEdge[j];
                    internal::addEdgesIntersectionToList(edge1IntersectionPoints, edge2IntersectionPoints, vertices,
                                                         vertexIndices[i], vertexIndices[(i + 1) % vertexIndices.size()],
                                                         vertexIndices[j], vertexIndices[(j + 1) % vertexIndices.size()]);
                }
            }
        }
        static void addPolygonsIntersectionsToList(std::vector<glm::dvec2>& vertices,
                                                   const std::vector<int>& polygon1VertexIndices,
                                                   std::vector<std::vector<internal::IntersectionPoint>>& intersectionPointsByEdgeForPolygon1,
                                                   const std::vector<int>& polygon2VertexIndices,
                                                   std::vector<std::vector<internal::IntersectionPoint>>& intersectionPointsByEdgeForPolygon2) {
            for (int i = 0; i < polygon1VertexIndices.size(); i++) {
                std::vector<internal::IntersectionPoint>& edge1IntersectionPoints = intersectionPointsByEdgeForPolygon1[i];
                for (int j = 0; j < polygon2VertexIndices.size(); j++) {
                    std::vector<internal::IntersectionPoint>& edge2IntersectionPoints = intersectionPointsByEdgeForPolygon2[j];
                    internal::addEdgesIntersectionToList(edge1IntersectionPoints, edge2IntersectionPoints, vertices,
                                                         polygon1VertexIndices[i], polygon1VertexIndices[(i + 1) % polygon1VertexIndices.size()],
                                                         polygon2VertexIndices[j], polygon2VertexIndices[(j + 1) % polygon2VertexIndices.size()]);
                }
            }
        }


        static void restoreLinkage(const std::vector<glm::dvec2>& vertices, const std::vector<int>& vertexIndices,
                                   const std::vector<std::vector<internal::IntersectionPoint>>& intersectionPointsByEdge,
                                   std::vector<VertexLinkage>& resultVertexLinkages) {
            int currentVertexIndex = -1;
            glm::dvec2 lastVertex = vertices[vertexIndices[0]], currentDirection(0, 0);
            for (int i = 0; i < vertexIndices.size(); i++) {
                if(currentVertexIndex != -1) resultVertexLinkages[currentVertexIndex].linkNext(vertexIndices[i], currentDirection);
                currentVertexIndex = vertexIndices[i];
                glm::dvec2 currentVertex = vertices[vertexIndices[(i + 1) % vertexIndices.size()]];
                currentDirection = glm::normalize(currentVertex - lastVertex);
                for(const internal::IntersectionPoint& intersectionPoint : intersectionPointsByEdge[i]) {
                    // We check this to exclude repeating intersectionPoints which may occur because of floating point errors
                    if(currentVertexIndex != intersectionPoint.vertexIndex) {
                        if(currentVertexIndex != -1) resultVertexLinkages[currentVertexIndex].linkNext(intersectionPoint.vertexIndex, currentDirection);
                        currentVertexIndex = intersectionPoint.vertexIndex;
                    }
                }
                lastVertex = currentVertex;
            }
            resultVertexLinkages[currentVertexIndex].linkNext(vertexIndices[0], currentDirection);
        }


        static std::vector<std::vector<std::vector<internal::IntersectionPoint>>> insertSteinerVertices(std::vector<glm::dvec2>& vertices,
                                                                                                        const std::vector<std::vector<int>>& polygonVertexIndices) {
            std::vector<std::vector<std::vector<internal::IntersectionPoint>>> intersectionPointsByEdge(polygonVertexIndices.size());
            for (int i = 0; i < polygonVertexIndices.size(); i++) {
                intersectionPointsByEdge[i] = std::vector<std::vector<internal::IntersectionPoint>>(polygonVertexIndices[i].size());
            }
            for (int i = 0; i < polygonVertexIndices.size(); i++) {
                internal::addPolygonSelfIntersectionsToList(vertices, polygonVertexIndices[i], intersectionPointsByEdge[i]);
            }
            for (int i = 0; i < (int) polygonVertexIndices.size() - 1; i++) {
                const std::vector<int>& polygon1VertexIndices = polygonVertexIndices[i];
                std::vector<std::vector<internal::IntersectionPoint>>& intersectionPointsByEdgeForPolygon1 = intersectionPointsByEdge[i];
                for (int j = i + 1; j < polygonVertexIndices.size(); j++) {
                    const std::vector<int>& polygon2VertexIndices = polygonVertexIndices[j];
                    std::vector<std::vector<internal::IntersectionPoint>>& intersectionPointsByEdgeForPolygon2 = intersectionPointsByEdge[j];
                    internal::addPolygonsIntersectionsToList(vertices,
                                                             polygon1VertexIndices, intersectionPointsByEdgeForPolygon1,
                                                             polygon2VertexIndices, intersectionPointsByEdgeForPolygon2);
                }
            }
            return intersectionPointsByEdge;
        }


    }


    DECOMPOSITION_API std::vector<VertexLinkage> insertSteinerVerticesForPolygons(std::vector<glm::dvec2>& vertices,
                                                                                  const std::vector<std::vector<int>>& polygonVertexIndices) {
        std::vector<std::vector<std::vector<internal::IntersectionPoint>>> intersectionPointsByEdge =
                internal::insertSteinerVertices(vertices, polygonVertexIndices);
        std::vector<VertexLinkage> resultVertexLinkages(vertices.size());
        for (int i = 0; i < polygonVertexIndices.size(); i++) {
            internal::restoreLinkage(vertices, polygonVertexIndices[i], intersectionPointsByEdge[i], resultVertexLinkages);
        }
        return resultVertexLinkages;
    }


}