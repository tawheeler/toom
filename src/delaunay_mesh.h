#ifndef DELAUNAY_MESH_H_INCLUDED
#define DELAUNAY_MESH_H_INCLUDED

#include "typedefs.h"
#include "vec.h"

#define INVALID_MESH_INDEX -1

// Determine whether the given points are in right-hand order (CCW).
// Returns a value > 0 if right-hand (CCW).
// Returns a value < 0 if left-hand (CW).
// Returns zero if the points are colinear.
f32 GetRightHandedness(const v2 *a, const v2 *b, const v2 *c);

// A quarter edge represents a directed edge in a QuadEdgeTree.camera->
// Each undirected edge in the graph A <-> B with face L on the left of A->B and face R on the right
// of A->B has four quarter edges:
//    A->B, L->R, B->A, and R->L.
// If this quarter edge is primal, then the vertex index points to its vertex.
// If this quarter edge is dual (i.e., represents a face), then its is set to typemax(size_t).
typedef struct QuarterEdge_s
{
    int index;
    v2 *vertex;                 // Non-null if this is a quarter edge originating at a vertex
    struct QuarterEdge_s *next; // The next right-hand (CCW) QuarterEdge with the same origin
    struct QuarterEdge_s *rot;  // The next right-hand (CCW) QuarterEdge associated with the same
                                // undirected original edge.
} QuarterEdge;

// A dual edge connects two faces.
bool IsDualEdge(const QuarterEdge *qe);

// A primal edge connects two vertices.
bool IsPrimalEdge(const QuarterEdge *qe);

// QuarterEdge traversal.
QuarterEdge *QENext(const QuarterEdge *qe);
QuarterEdge *QERot(const QuarterEdge *qe);
QuarterEdge *QESym(const QuarterEdge *qe);
QuarterEdge *QETor(const QuarterEdge *qe);
QuarterEdge *QEPrev(const QuarterEdge *qe);
QuarterEdge *QELnext(const QuarterEdge *qe);

// A DelaunayMesh is a Quad Edge Mesh that can dynamically generate a Delaunay (or constrained
// Delaunay) triangulation. A Quad Edge Mesh's name comes from the fact that each edge in a "normal"
// mesh graph that connects two vertices is represented by 4 edges in the quad edge tree - one
// directed edge between each vertex and one directed edge between the two adjoined faces. Because
// of this, we always have a multiple of 4 quater edges.
//
// This class allows for the insertion and removal of items. As such, our vectors may change in
// length as memory reallocates. We use integers to index into our vertices and quarter edges to
// avoid memeory reallocation issues.
struct DelaunayMesh
{
    // The square of the maximum radius that a point can be from the origin
    f32 square_bounding_radius;

    // Points within this distance are considered colinear
    f32 square_min_dist_to_vertex;

    // A point this close to an edge is considered coincident
    f32 min_dist_to_edge;

    // All of the stored vertices.
    // This class owns this memory.
    u32 n_vertices;
    u32 max_n_vertices;
    v2 *vertices;

    // All of the stored quarter edges.
    // This class owns this memory.
    u32 n_quarter_edges;
    u32 max_n_quarter_edges;
    QuarterEdge *quarter_edges;
};

// Allocate a DelaunayMesh and initialize it to be a bounding triangle with
// the given radius.
struct DelaunayMesh *ConstructEmptyDelaunayMesh(
    f32 bounding_radius,
    f32 min_dist_to_vertex,
    f32 min_dist_to_edge,
    u32 max_n_vertices,
    u32 max_n_quarter_edges);

// Deallocate and free the given mesh.
void DeconstructDelaunayMesh(struct DelaunayMesh *mesh);

size_t DelaunayMeshNumVertices(struct DelaunayMesh *mesh);
size_t DelaunayMeshNumQuarterEdges(struct DelaunayMesh *mesh);
size_t DelaunayMeshNumEdges(struct DelaunayMesh *mesh);

const v2 *DelaunayMeshGetVertex(struct DelaunayMesh *mesh, int i);
QuarterEdge *DelaunayMeshGetQuarterEdge(struct DelaunayMesh *mesh, int i);

// Whether the given vertex is one of the boundary vertices (The first three in vertices_).
bool DelaunayMeshIsBoundaryVertex(const struct DelaunayMesh *mesh, const v2 *v);

// Get the quarter edge pointing from i to j, if it exists.
// Note that the current implementation loops over all quarter edges.
QuarterEdge *DelaunayMeshGetQuarterEdgeBetween(struct DelaunayMesh *mesh, int i, int j);

// Returns true if an edge from i to j exists.
// Note that the current implementation loops over all quarter edges.
bool DelaunayMeshHasEdge(struct DelaunayMesh *mesh, int i, int j);

// Insert a new vertex into the mesh and it to continue to be a Delaunay triangularization.
// Returns the index of the added vertex if a new point was added, and kInvalidIndex instead.
// (If the new point is coincident with an existing point or the boundary edge, no new point is added.)
int DelaunayMeshAddVertex(struct DelaunayMesh *mesh, const v2 *p);

// Enforce an edge between the ith and jth vertices.
// Returns true if this operation was successful.
bool DelaunayMeshConstrainEdge(struct DelaunayMesh *mesh, int i, int j);

// Find the triangle in our mesh that encloses the given point.
// The search starts from the given dual quarter edge.
// Return a dual quarter edge originating from the triangle face.
QuarterEdge *DelaunayMeshGetEnclosingTriangle(struct DelaunayMesh *mesh, const v2 *p, QuarterEdge *qe_dual);
QuarterEdge *DelaunayMeshGetEnclosingTriangle2(struct DelaunayMesh *mesh, const v2 *p);

// Given a dual quarter edge, return the first vertex on that face.
const v2 *DelaunayMeshGetTriangleVertex1(struct DelaunayMesh *mesh, const QuarterEdge *qe_dual);
const v2 *DelaunayMeshGetTriangleVertex2(struct DelaunayMesh *mesh, const QuarterEdge *qe_dual);
const v2 *DelaunayMeshGetTriangleVertex3(struct DelaunayMesh *mesh, const QuarterEdge *qe_dual);

#endif