#include "gtest/gtest.h"
#define DECOMPOSITION_TEST

#include "edge-intersection.cpp"

using namespace decomposition;
using namespace decomposition::internal;



static std::vector<VertexLinkage> buildPolygonGraph(std::vector<glm::dvec2>& vertices) {
    std::vector<int> indices(vertices.size());
    for (int i = 0; i < vertices.size(); i++) indices[i] = i;
    return insertSteinerVerticesForPolygons(vertices, {indices});
}



TEST(EdgeIntersection, BinarySearch) {
    std::vector<IntersectionPoint> points {
            {0, 0.0},
            {0, 0.1},
            {0, 0.5},
            {0, 0.9},
            {0, 1.0}
    };
    ASSERT_EQ(binarySearchInIntersectionPoints(points, 0), 0);
    ASSERT_EQ(binarySearchInIntersectionPoints(points, 0.05), 1);
    ASSERT_EQ(binarySearchInIntersectionPoints(points, 0.1), 1);
    ASSERT_EQ(binarySearchInIntersectionPoints(points, 0.2), 2);
    ASSERT_EQ(binarySearchInIntersectionPoints(points, 0.3), 2);
    ASSERT_EQ(binarySearchInIntersectionPoints(points, 0.4), 2);
    ASSERT_EQ(binarySearchInIntersectionPoints(points, 0.5), 2);
    ASSERT_EQ(binarySearchInIntersectionPoints(points, 0.6), 3);
    ASSERT_EQ(binarySearchInIntersectionPoints(points, 0.7), 3);
    ASSERT_EQ(binarySearchInIntersectionPoints(points, 0.8), 3);
    ASSERT_EQ(binarySearchInIntersectionPoints(points, 0.9), 3);
    ASSERT_EQ(binarySearchInIntersectionPoints(points, 0.95), 4);
    ASSERT_EQ(binarySearchInIntersectionPoints(points, 1.0), 4);
}


TEST(EdgeIntersection, SimpleIntersectingEdges) {
    {
        auto intersection = findIntersectionPoint(glm::dvec2(-1, -1), glm::dvec2(1, 1), glm::dvec2(-1, 1), glm::dvec2(1, -1));
        ASSERT_TRUE(!!intersection);
        ASSERT_DOUBLE_EQ(intersection->point.x, 0);
        ASSERT_DOUBLE_EQ(intersection->point.y, 0);
        ASSERT_DOUBLE_EQ(intersection->edge1Offset, 0.5);
        ASSERT_DOUBLE_EQ(intersection->edge2Offset, 0.5);
    }
    {
        auto intersection = findIntersectionPoint(glm::dvec2(-1, 0), glm::dvec2(3, 0), glm::dvec2(2, 3), glm::dvec2(2, -1));
        ASSERT_TRUE(!!intersection);
        ASSERT_DOUBLE_EQ(intersection->point.x, 2);
        ASSERT_DOUBLE_EQ(intersection->point.y, 0);
        ASSERT_DOUBLE_EQ(intersection->edge1Offset, 0.75);
        ASSERT_DOUBLE_EQ(intersection->edge2Offset, 0.75);
    }
}


TEST(EdgeIntersection, TouchingEdges) {
    // We must treat segments as [startVertex, endVertex), so collisions between edge ends must not be detected
    ASSERT_TRUE (!!findIntersectionPoint(glm::dvec2(7654, 213), glm::dvec2(325, 20), glm::dvec2(7654, 213), glm::dvec2(325, 0)));
    ASSERT_TRUE (!!findIntersectionPoint(glm::dvec2(7654, 213), glm::dvec2(325, 0), glm::dvec2(7654, 213), glm::dvec2(325, 20)));
    ASSERT_FALSE(!!findIntersectionPoint(glm::dvec2(7654, 213), glm::dvec2(325, 20), glm::dvec2(325, 0), glm::dvec2(7654, 213)));
    ASSERT_FALSE(!!findIntersectionPoint(glm::dvec2(325, 0), glm::dvec2(7654, 213), glm::dvec2(7654, 213), glm::dvec2(325, 20)));

    ASSERT_TRUE (!!findIntersectionPoint(glm::dvec2(0, -1), glm::dvec2(0, 1), glm::dvec2(0, 0), glm::dvec2(1, 0)));
    ASSERT_FALSE(!!findIntersectionPoint(glm::dvec2(0, -1), glm::dvec2(0, 1), glm::dvec2(0, 1), glm::dvec2(0, 0)));

    ASSERT_TRUE (!!findIntersectionPoint(glm::dvec2(0, 0), glm::dvec2(1, 0), glm::dvec2(0, -1), glm::dvec2(0, 1)));
    ASSERT_FALSE(!!findIntersectionPoint(glm::dvec2(0, 1), glm::dvec2(0, 0), glm::dvec2(0, -1), glm::dvec2(0, 1)));
}


TEST(EdgeIntersection, OnePointIntersecionOf3Edges) {
    std::vector<glm::dvec2> vertices {
        {-1, 0},
        {1, 0},
        {1, 1},
        {-1, -1},
        {0, -1},
        {0, 1}
    };
    auto linkage = buildPolygonGraph(vertices);
    ASSERT_EQ(vertices.size(), 7);
    ASSERT_DOUBLE_EQ(vertices[6].x, 0);
    ASSERT_DOUBLE_EQ(vertices[6].y, 0);
    ASSERT_EQ(linkage[6].pathsNumber, 3);
}


TEST(EdgeIntersection, OnePointIntersecionOf4Edges) {
    std::vector<glm::dvec2> vertices {
            {-10, 0},
            {1, 0},
            {1, 1},
            {-1, -1},
            {0, -1},
            {0, 1},
            {-1, 1},
            {10, -10}
    };
    auto linkage = buildPolygonGraph(vertices);
    ASSERT_EQ(vertices.size(), 9);
    ASSERT_DOUBLE_EQ(vertices[8].x, 0);
    ASSERT_DOUBLE_EQ(vertices[8].y, 0);
    ASSERT_EQ(linkage[8].pathsNumber, 4);
}


TEST(EdgeIntersection, SelfTouchingPolygon) {
    // We must not create additional vertices which duplicate existing ones
    std::vector<glm::dvec2> vertices {
            {-1, 0},
            {1, 0},
            {1, 1},
            {0, 0},
            {-1, 1}
    };
    auto linkage = buildPolygonGraph(vertices);
    ASSERT_EQ(vertices.size(), 5);
    ASSERT_EQ(linkage[3].pathsNumber, 2);
}