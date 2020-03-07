#include "gtest/gtest.h"
#define DECOMPOSITION_TEST

#include "edge-intersection.cpp"
#include "polygon-decomposition.cpp"

using namespace decomposition;
using namespace decomposition::internal;



static std::vector<Polygon> decomposePolygon(std::vector<glm::dvec2>& vertices) {
    std::vector<int> indices(vertices.size());
    for (int i = 0; i < vertices.size(); i++) indices[i] = i;
    return decomposePolygonGraph(vertices, insertSteinerVerticesForPolygons(vertices, {indices}));
}



TEST(PolygonDecomposition, PathChoosing) {
    for (int order = 0; order < 2; order++) {
        VertexLinkage vertex;
        int index = 0;
        for (double angle = 0.1; angle <= 6.2; angle += 0.1) {
            vertex.linkNext(index++,  glm::dvec2(-sin(angle), -cos(angle)));
        }
        for (int i = 0; i < index; i++) {
            VertexLinkage::Path& path = choosePath(vertex, glm::dvec2(0, 1),
                    order == 0 ? WindingOrder::COUNTER_CLOCKWISE : WindingOrder::CLOCKWISE);
            ASSERT_EQ(path.next, order == 0 ? i : index - i - 1);
            path.next = -1;
        }
    }
}


TEST(PolygonDecomposition, SimpleSelfIntersection) {
    std::vector<glm::dvec2> vertices {
            {-1, -1},
            {-1, 1},
            {1, -1},
            {1, 1}
    };
    auto polygons = decomposePolygon(vertices);
    ASSERT_EQ(polygons.size(), 2);
    ASSERT_EQ(polygons[0].getWinding(), -polygons[1].getWinding());
}


TEST(PolygonDecomposition, ComplexSelfIntersection) {
    std::vector<glm::dvec2> vertices {
            {-30, 20},
            {0, -30},
            {0, 1},
            {1, 1},
            {-1, -1},
            {30, -1},
            {0, 30},
            {-30, -20}
    };
    auto polygons = decomposePolygon(vertices);
    ASSERT_EQ(vertices.size(), 11);
    ASSERT_EQ(polygons.size(), 4);
    int windingsSum = 0;
    for (Polygon& polygon : polygons) windingsSum += polygon.getWinding();
    ASSERT_EQ(windingsSum, 0); // There must be 2 clockwise and 2 counter-clockwise polygons
}