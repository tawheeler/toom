#ifndef DELAUNAY_MESH_H_INCLUDED
#define DELAUNAY_MESH_H_INCLUDED

#include "typedefs.h"
#include "vec.h"

#define INVALID_MESH_INDEX -1

// A quarter edge represents a directed edge in a QuadEdgeTree.
// Each undirected edge in the graph A <-> B with face L on the left of A->B and face R on the right
// of A->B has four quarter edges:
//    A->B, L->R, B->A, and R->L.
// If this quarter edge is primal, then the vertex index points to its vertex.
// If this quarter edge is dual (i.e., represents a face), then its is set to typemax(size_t).
struct QuarterEdge {
    v2* vertex;        // Non-null if this is a quarter edge originating at a vertex
    struct QuarterEdge* next; // The next right-hand (CCW) QuarterEdge with the same origin
    struct QuarterEdge* rot;  // The next right-hand (CCW) QuarterEdge associated with the same
                       // undirected original edge.
};

// A dual edge connects two faces.
bool IsDualEdge(const struct QuarterEdge* qe);

// A primal edge connects two vertices.
bool IsPrimalEdge(const struct QuarterEdge* qe);

// QuarterEdge traversal.
struct QuarterEdge* QENext(const struct QuarterEdge* qe);
struct QuarterEdge* QERot(const struct QuarterEdge* qe);
struct QuarterEdge* QESym(const struct QuarterEdge* qe);
struct QuarterEdge* QETor(const struct QuarterEdge* qe);
struct QuarterEdge* QEPrev(const struct QuarterEdge* qe);
struct QuarterEdge* QELnext(const struct QuarterEdge* qe);

// A DelaunayMesh is a Quad Edge Mesh that can dynamically generate a Delaunay (or constrained
// Delaunay) triangulation. A Quad Edge Mesh's name comes from the fact that each edge in a "normal"
// mesh graph that connects two vertices is represented by 4 edges in the quad edge tree - one
// directed edge between each vertex and one directed edge between the two adjoined faces. Because
// of this, we always have a multiple of 4 quater edges.
//
// This class allows for the insertion and removal of items. As such, our vectors may change in
// length as memory reallocates. We use integers to index into our vertices and quarter edges to
// avoid memeory reallocation issues.
struct DelaunayMesh {
    // The maximum radius that a point can be from the origin
    f32 bounding_radius;

    // Points within this distance are considered colinear
    f32 min_dist_to_vertex;

    // A point this close to an edge is considered coincident
    f32 min_dist_to_edge;

    // All of the stored vertices.
    // This class owns this memory.
    u32 n_vertices;
    u32 max_n_vertices;
    v2* vertices;

    // All of the stored quarter edges.
    // This class owns this memory.
    u32 n_quarter_edges;
    u32 max_n_quarter_edges;
    struct QuarterEdge* quarter_edges;
};

// Allocate a DelaunayMesh and initialize it to be a bounding triangle with
// the given radius.
struct DelaunayMesh* ConstructEmptyDelaunayMesh(
    f32 bounding_radius,
    f32 min_dist_to_vertex,
    f32 min_dist_to_edge,
    u32 max_n_vertices,
    u32 max_n_quarter_edges);

// Deallocate and free the given mesh. 
void DeconstructDelaunayMesh(struct DelaunayMesh* mesh);

size_t DelaunayMeshNumVertices(struct DelaunayMesh* mesh);
size_t DelaunayMeshNumQuarterEdges(struct DelaunayMesh* mesh);
size_t DelaunayMeshNumEdges(struct DelaunayMesh* mesh);

const v2* DelaunayMeshGetVertex(struct DelaunayMesh* mesh, int i);
struct QuarterEdge* DelaunayMeshGetQuarterEdge(struct DelaunayMesh* mesh, int i);

// Get the quarter edge pointing from i to j, if it exists.
// Note that the current implementation loops over all quarter edges.
struct QuarterEdge* DelaunayMeshGetQuarterEdgeBetween(struct DelaunayMesh* mesh, int i, int j);

// Returns true if an edge from i to j exists.
// Note that the current implementation loops over all quarter edges.
bool DelaunayMeshHasEdge(struct DelaunayMesh* mesh, int i, int j);

// Add a new vertex to the mesh, returning the index of the newly added vertex.
// Returns NULL if we are at capacity.
v2* _DelaunayMeshAddVertex(struct DelaunayMesh* mesh, f32 x, f32 y);

// Add a new edge between the vertices at index a and b.
// Create the quarter edges associated with the given undirected edge.
// We always create quarter-edges in groups of four.
// Returns a reference to the quarter-edge from A to B.
// Returns NULL if we are at capacity.
struct QuarterEdge* _DelaunayMeshAddEdge(struct DelaunayMesh* mesh, v2* a, v2* b);

// A utility function used to join quarter edges.
void _DelaunayMeshSplice(struct DelaunayMesh* mesh, struct QuarterEdge* a, struct QuarterEdge* b);

#endif