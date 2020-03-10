# decomposition-library
##### Small and lightweight C++ library for polygon clipping and triangulation.

It's split into few independent parts, allowing for easy testing and upgrading some pieces without affecting other ones.

###### Features:
1. Finding intersections of polygons (including self-intersections)
2. Decomposing polygons into set of simple, non-intersecting shapes
3. Reconstructing hierarchy of shapes and polygons with holes
4. Triangulation of simple polygons with holes

## How to install
1. `git clone https://github.com/YaaZ/decomposition-library.git`
2. `git submodule update --init --recursive`
3. Add `add_subdirectory(<path-to-library-directory>)` in your project  Cmake
4. Add `target_link_libraries(<your-project-target> decomposition_library)` in your project  Cmake

## How to use
API of this library is split into few functions.
Almost all functions consume data, produced by previous steps for performance reasons, in most cases
you will need to do copy of data yourself if you want to use it somewhere else.
See more details about functions usage in [doxygen comments](include/decomposition.h).

#### 1. Inserting Steiner vertices
```cpp
decomposition::insertSteinerVerticesForPolygons(vertexCoordinates, vertexIndices);
```
Will find all intersections of polygons (including self-intersections), insert them into polygons as
additional vertices and build graph representing new polygons.

#### 2. Decomposing polygon graph
```cpp
decomposition::decomposePolygonGraph(vertexCoordinates, polygonGraph)
```
Will consume graph, created on previous step and decompose it into set of simple, non-intersecting polygons.
It uses winding order of polygons to determine, wheter some polygon overlaps another polygon, or it represents a hole.
In short: polygons with same winding order sums up and polygon with opposite winding order subtracts.
For better understanding you can try [decomposition-viewer utility](https://github.com/YaaZ/decomposition-viewer).

#### 3. Building polygon trees
```cpp
decomposition::buildPolygonTrees(vertexCoordinates, polygons)
```
Will consume polygons, generated on previous step and compose them into trees, representing hierarchy of nested shapes.
These trees can be used for purposes like fast point lookups, to determine, how many layers of polygons are under given point.

#### 4. Building polygon area trees
```cpp
decomposition::buildPolygonAreaTrees(polygonTrees)
```
Will consume polygon trees, built on previous step and rebuild them into trees, representing actual hierarchy of polygon areas
rather than hierarchy of its contours. That is, if on previous step there were trees of polygons that had to be summed
or subtracted based on their winding order, now there are only nested into each other polygons with holes, winding order
of which don't make as much sense as before.

#### 5. Triangulating polygons
```cpp
decomposition::triangulatePolygonWithHoles(vertexCoordinates, polygonWithHoles)
```
Will generate triangles from one of polygons with holes, created on previous step.

### Usage example
```cpp
std::vector<glm::dvec2> vertices {
        {0, 0},
        {1, 0},
        {0, 1}
};
const std::vector<std::vector<int>> vertexIndices {
        {0, 1, 2}
};

std::vector<decomposition::VertexLinkage> polygonVertexGraph =
        decomposition::insertSteinerVerticesForPolygons(vertices, vertexIndices);

std::vector<decomposition::Polygon> polygonSet =
        decomposition::decomposePolygonGraph(vertices, std::move(polygonVertexGraph));

std::vector<decomposition::PolygonTree> polygonTree =
        decomposition::buildPolygonTrees(vertices, std::move(polygonSet));

std::vector<decomposition::PolygonWithHolesTree> polygonAreaTree =
        decomposition::buildPolygonAreaTrees(std::move(polygonTree));

// Iterate roots only (do not triangulate overlapping areas more than once)
for(const decomposition::PolygonWithHoles& polygon : polygonAreaTree) {
    std::vector<glm::ivec3> polygonTriangles =
            decomposition::triangulatePolygonWithHoles(vertices, polygon);
}
```

# 

This library uses [GLM](https://github.com/g-truc/glm) and [Google Test](https://github.com/google/googletest) libraries
