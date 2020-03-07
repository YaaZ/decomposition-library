#pragma once


#include <cstdint>
#include <vector>
#include <glm.hpp>


namespace decomposition {


    enum WindingOrder {
        CLOCKWISE = 1,
        COUNTER_CLOCKWISE = -1
    };


    /**
     * Linkage object contains all edge paths intersecting some vertex. More specifically each linkage object contains
     * set of paths: vertex indices, so that there is an edge from current vertex to vertex denoted by each of its paths
     */
    class VertexLinkage {
    public:
        struct Path {
            int next;
            /// Normalized vector pointing to next vertex
            glm::dvec2 directionToNextVertex;
        };

    private:
        /*
         * We split paths into 2 pieces: static array and vector, we do this to avoid frequent dynamic allocations.
         * Many intersections in one point are not likely, so most of the time there will be no more that 2 paths.
         */
        static constexpr size_t PATHS_STATIC_PART_SIZE = 2;
        Path pathsStaticPart[PATHS_STATIC_PART_SIZE];
        std::vector<Path> pathsDynamicPart;
    public:
        int pathsNumber {0};
        inline VertexLinkage() {};
        inline Path& getPath(int i) {
            if(i < PATHS_STATIC_PART_SIZE) return pathsStaticPart[i];
            else return pathsDynamicPart[i - PATHS_STATIC_PART_SIZE];
        }
        inline void linkNext(int nextVertexIndex, glm::dvec2 edgeDirection) {
            if(pathsNumber < PATHS_STATIC_PART_SIZE) pathsStaticPart[pathsNumber] = {nextVertexIndex, edgeDirection};
            else pathsDynamicPart.push_back({nextVertexIndex, edgeDirection});
            pathsNumber++;
        }

    };


    struct Polygon {

        std::vector<int> vertexIndices;
        struct AABB {
            glm::dvec2 min, max;
            inline void expand(glm::dvec2 point) {
                if(point.x < min.x) min.x = point.x;
                if(point.y < min.y) min.y = point.y;
                if(point.x > max.x) max.x = point.x;
                if(point.y > max.y) max.y = point.y;
            }
        } aabb;
        double signedArea;

        inline double getArea() const { return glm::abs(signedArea); }
        inline WindingOrder getWinding() const { return signedArea >= 0 ? WindingOrder::CLOCKWISE : WindingOrder::COUNTER_CLOCKWISE; }

        inline bool operator<(const Polygon& p) { return getArea() < p.getArea(); }
        inline bool operator<=(const Polygon& p) { return getArea() <= p.getArea(); }
        inline bool operator>(const Polygon& p) { return getArea() > p.getArea(); }
        inline bool operator>=(const Polygon& p) { return getArea() >= p.getArea(); }

    };


    struct PolygonTree : public Polygon {

        std::vector<PolygonTree> childrenPolygons;
        /**
         * Total winding: sum of winding values of current vertex and all its parents.
         * Its absolute value has a meaning of number of layers, or number of polygons intersecting in this area.
         * Therefore net winding of 0 means that this is a "hole" polygon
         */
        int netWinding;

    };


    struct PolygonWithHoles : public Polygon {

        /**
         * "Hole" polygons, must have opposite winding to winding of current polygon
         */
        std::vector<Polygon> holes;

    };

    struct PolygonWithHolesTree : public PolygonWithHoles {

        std::vector<PolygonWithHolesTree> childrenPolygons;
        /// @see PolygonTree.netWinding
        int netWinding;

    };


}