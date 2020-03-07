#include "gtest/gtest.h"
#define DECOMPOSITION_TEST

#include "edge-intersection.cpp"
#include "polygon-decomposition.cpp"
#include "polygon-tree.cpp"

using namespace decomposition;
using namespace decomposition::internal;



static std::vector<PolygonTree> createPolygonTree(std::vector<glm::dvec2>& vertices) {
    std::vector<int> indices(vertices.size());
    for (int i = 0; i < vertices.size(); i++) indices[i] = i;
    return buildPolygonTrees(vertices, decomposePolygonGraph(vertices, insertSteinerVerticesForPolygons(vertices, {indices})));
}
static PolygonTree createPolygonTree(WindingOrder winding, int netWinding) {
    PolygonTree tree;
    tree.signedArea = (double) winding;
    tree.netWinding = netWinding;
    return tree;
}



TEST(PolygonTree, PointToPolygonNesting) {
    std::vector<glm::dvec2> vertices {
            {-1, -1},
            {-1, 1},
            {0, 0},
            {1, 1},
            {1, -1}
    };

    std::vector<int> indices[2] {{0, 1, 2, 3, 4}, {4, 3, 2, 1, 0}}; // We check both directions for the polygon

    for (int i = 0; i < 2; i++) {
        ASSERT_EQ(calculatePointToPolygonNesting(vertices, indices[i], glm::dvec2(-0.75, 0.5)), 1);
        ASSERT_EQ(calculatePointToPolygonNesting(vertices, indices[i], glm::dvec2(0.75, 0.5)), 1);
        ASSERT_EQ(calculatePointToPolygonNesting(vertices, indices[i], glm::dvec2(-0.5, 0)), 1);
        ASSERT_EQ(calculatePointToPolygonNesting(vertices, indices[i], glm::dvec2(0.5, 0)), 1);

        ASSERT_EQ(calculatePointToPolygonNesting(vertices, indices[i], glm::dvec2(-1, 0)), 0);
        ASSERT_EQ(calculatePointToPolygonNesting(vertices, indices[i], glm::dvec2(1, 0)), 0);
        ASSERT_EQ(calculatePointToPolygonNesting(vertices, indices[i], glm::dvec2(-0.5, 0.5)), 0);
        ASSERT_EQ(calculatePointToPolygonNesting(vertices, indices[i], glm::dvec2(0.5, 0.5)), 0);
        ASSERT_EQ(calculatePointToPolygonNesting(vertices, indices[i], glm::dvec2(0, 0)), 0);
        ASSERT_EQ(calculatePointToPolygonNesting(vertices, indices[i], glm::dvec2(0, -1)), 0);
        ASSERT_EQ(calculatePointToPolygonNesting(vertices, indices[i], glm::dvec2(-1, -1)), 0);
        ASSERT_EQ(calculatePointToPolygonNesting(vertices, indices[i], glm::dvec2(-1, 1)), 0);
        ASSERT_EQ(calculatePointToPolygonNesting(vertices, indices[i], glm::dvec2(1, 1)), 0);
        ASSERT_EQ(calculatePointToPolygonNesting(vertices, indices[i], glm::dvec2(1, -1)), 0);

        ASSERT_EQ(calculatePointToPolygonNesting(vertices, indices[i], glm::dvec2(0, 0.5)), -1);
        ASSERT_EQ(calculatePointToPolygonNesting(vertices, indices[i], glm::dvec2(0, 2)), -1);
        ASSERT_EQ(calculatePointToPolygonNesting(vertices, indices[i], glm::dvec2(-2, 0)), -1);
        ASSERT_EQ(calculatePointToPolygonNesting(vertices, indices[i], glm::dvec2(2, 0)), -1);
    }
}


TEST(PolygonTree, ComplexSelfIntersection) {
    std::vector<glm::dvec2> vertices {
            {-2, -2},
            {-2, 2},
            {10, 2},
            {-1, 1},
            {1, -1},
            {2, 10},
            {2, -2}
    };
    auto tree = createPolygonTree(vertices);
    ASSERT_EQ(tree.size(), 1);
    ASSERT_EQ(tree[0].childrenPolygons.size(), 2);
    ASSERT_EQ(tree[0].netWinding, 1);
    ASSERT_EQ(tree[0].childrenPolygons[0].childrenPolygons.size(), 0);
    ASSERT_EQ(tree[0].childrenPolygons[1].childrenPolygons.size(), 0);
    ASSERT_EQ(tree[0].childrenPolygons[0].netWinding + tree[0].childrenPolygons[1].netWinding, 2);
}


TEST(PolygonTree, RealCases) { // Edge cases which were known to have bugs or cause crashes
    {
        std::vector<glm::dvec2> vertices {
                {-1, -1},
                {1, -1},
                {1, 1},
                {-1, 1},
                {0, -0.75},
                {0, 0.5},
                {1, 0}
        };
        auto tree = createPolygonTree(vertices);
        ASSERT_EQ(tree.size(), 1);
        ASSERT_EQ(tree[0].netWinding, -1);
        ASSERT_EQ(tree[0].childrenPolygons.size(), 2);
        ASSERT_EQ(tree[0].childrenPolygons[0].childrenPolygons.size(), 0);
        ASSERT_EQ(tree[0].childrenPolygons[1].childrenPolygons.size(), 0);
        ASSERT_EQ(tree[0].childrenPolygons[0].netWinding + tree[0].childrenPolygons[1].netWinding, -2);
    }
    {
        std::vector<glm::dvec2> vertices {
                {441, 174},
                {296, 325},
                {440, 333},
                {403, 170},
                {270, 296},
                {236, 106},
                {580, 66},
                {619, 352},
                {245, 401},
                {193, 89},
                {474, 55},
                {540, 414},
                {666, 432},
                {662, 16},
                {215, 40},
                {264, 303},
                {470, 311}
        };
        auto tree = createPolygonTree(vertices);
        ASSERT_EQ(tree.size(), 3);
    }
    {
        std::vector<glm::dvec2> vertices {
                {560, 128},
                {342, 117},
                {560, 243},
                {499, 156},
                {560, 360}
        };
        auto tree = createPolygonTree(vertices);
        ASSERT_EQ(tree.size(), 1);
    }
    {
        std::vector<glm::dvec2> vertices {
                {239, 117},
                {186, 377},
                {394, 428},
                {629, 357},
                {513, 95},
                {481, 137},
                {553, 335},
                {389, 378},
                {238, 346},
                {281, 153},
                {277, 344},
                {286, 186},
                {458, 151},
                {515, 337},
                {257, 312},
                {315, 161},
                {479, 163},
                {539, 323},
                {365, 359},
                {257, 333},
                {296, 167},
                {545, 244},
                {416, 354}
        };
        auto tree = buildPolygonTrees(vertices, decomposePolygonGraph(vertices, insertSteinerVerticesForPolygons(vertices, {
                {0, 1, 2, 3, 4},
                {5, 6, 7, 8, 9},
                {10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21}
        })));
        ASSERT_EQ(tree.size(), 1);
        ASSERT_EQ(tree[0].childrenPolygons.size(), 1);
        ASSERT_EQ(tree[0].childrenPolygons[0].childrenPolygons.size(), 1);
        ASSERT_EQ(tree[0].childrenPolygons[0].childrenPolygons[0].childrenPolygons.size(), 1);
        ASSERT_EQ(tree[0].childrenPolygons[0].childrenPolygons[0].childrenPolygons[0].childrenPolygons.size(), 1);
        ASSERT_EQ(tree[0].childrenPolygons[0].childrenPolygons[0].childrenPolygons[0].childrenPolygons[0].childrenPolygons.size(), 0);
    }
}


TEST(PolygonTree, PolygonWithHolesTreeBuilding) {
    std::vector<PolygonTree> trees {createPolygonTree(CLOCKWISE, 1)};
    trees[0].childrenPolygons.push_back(createPolygonTree(CLOCKWISE, 2));
    trees[0].childrenPolygons[0].childrenPolygons.push_back(createPolygonTree(COUNTER_CLOCKWISE, 1));
    trees[0].childrenPolygons[0].childrenPolygons[0].childrenPolygons.push_back(createPolygonTree(CLOCKWISE, 2));
    trees[0].childrenPolygons[0].childrenPolygons[0].childrenPolygons.push_back(createPolygonTree(COUNTER_CLOCKWISE, 0));
    trees[0].childrenPolygons[0].childrenPolygons[0].childrenPolygons[1].childrenPolygons.push_back(createPolygonTree(COUNTER_CLOCKWISE, -1));
    trees[0].childrenPolygons[0].childrenPolygons.push_back(createPolygonTree(CLOCKWISE, 3));
    trees[0].childrenPolygons[0].childrenPolygons.push_back(createPolygonTree(CLOCKWISE, 3));
    trees[0].childrenPolygons[0].childrenPolygons[2].childrenPolygons.push_back(createPolygonTree(COUNTER_CLOCKWISE, 2));
    trees[0].childrenPolygons[0].childrenPolygons.push_back(createPolygonTree(COUNTER_CLOCKWISE, 1));
    trees[0].childrenPolygons[0].childrenPolygons[3].childrenPolygons.push_back(createPolygonTree(COUNTER_CLOCKWISE, 0));
    std::vector<PolygonWithHolesTree> result = buildPolygonAreaTrees(std::move(trees));
    ASSERT_EQ(result.size(), 2);
    ASSERT_EQ(result[0].netWinding, 1);
    ASSERT_EQ(result[0].holes.size(), 2);
    ASSERT_EQ(result[0].childrenPolygons.size(), 2);
    ASSERT_EQ(result[0].childrenPolygons[0].holes.size(), 2);
    ASSERT_EQ(result[0].childrenPolygons[0].childrenPolygons.size(), 2);
    ASSERT_EQ(result[0].childrenPolygons[0].childrenPolygons[0].holes.size(), 1);
    ASSERT_EQ(result[0].childrenPolygons[0].childrenPolygons[0].childrenPolygons.size(), 0);
    ASSERT_EQ(result[0].childrenPolygons[0].childrenPolygons[1].holes.size(), 0);
    ASSERT_EQ(result[0].childrenPolygons[0].childrenPolygons[1].childrenPolygons.size(), 0);
    ASSERT_EQ(result[0].childrenPolygons[1].holes.size(), 0);
    ASSERT_EQ(result[0].childrenPolygons[1].childrenPolygons.size(), 0);
    ASSERT_EQ(result[1].netWinding, -1);
    ASSERT_EQ(result[1].holes.size(), 0);
    ASSERT_EQ(result[1].childrenPolygons.size(), 0);
}