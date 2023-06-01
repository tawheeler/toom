#include "delaunay_mesh.h"

#include "constants.h"

#include <math.h>
#include <stdlib.h>

// ------------------------------------------------------------------------------------------------
bool IsDualEdge(const struct QuarterEdge* qe) { return qe->vertex == NULL; }

// ------------------------------------------------------------------------------------------------
bool IsPrimalEdge(const struct QuarterEdge* qe) { return qe->vertex != NULL; }

// ------------------------------------------------------------------------------------------------
struct QuarterEdge* QENext(const struct QuarterEdge* qe) { return qe->next; }
struct QuarterEdge* QERot(const struct QuarterEdge* qe) { return qe->rot; }
struct QuarterEdge* QESym(const struct QuarterEdge* qe) { return qe->rot->rot; }
struct QuarterEdge* QETor(const struct QuarterEdge* qe) { return qe->rot->rot->rot; }
struct QuarterEdge* QEPrev(const struct QuarterEdge* qe) { return qe->rot->next->rot; }
struct QuarterEdge* QELnext(const struct QuarterEdge* qe) { return QETor(qe)->next->rot; }

// ------------------------------------------------------------------------------------------------
void QESwapNexts(struct QuarterEdge* a, struct QuarterEdge* b) {
    struct QuarterEdge* tmp = a->next;
    a->next = b->next;
    b->next = tmp;
}

// ------------------------------------------------------------------------------------------------
void QESplice(struct QuarterEdge* a, struct QuarterEdge* b) {
    QESwapNexts(a->next->rot, b->next->rot);
    QESwapNexts(a, b);
}

// ------------------------------------------------------------------------------------------------
struct DelaunayMesh* ConstructEmptyDelaunayMesh(
    f32 bounding_radius,
    f32 min_dist_to_vertex,
    f32 min_dist_to_edge,
    u32 max_n_vertices,
    u32 max_n_quarter_edges) {
    
    struct DelaunayMesh* mesh = (struct DelaunayMesh*) malloc(sizeof(struct DelaunayMesh));
    mesh->bounding_radius = bounding_radius;
    mesh->min_dist_to_vertex = min_dist_to_vertex;
    mesh->min_dist_to_edge = min_dist_to_edge;

    mesh->max_n_quarter_edges = max_n_quarter_edges;
    mesh->max_n_vertices = max_n_vertices;
    mesh->n_quarter_edges = 0;
    mesh->n_vertices = 0;

    mesh->vertices = (v2*) malloc(max_n_vertices * sizeof(v2));
    mesh->quarter_edges = (struct QuarterEdge*) malloc(max_n_quarter_edges * sizeof(struct QuarterEdge));

    // The triangle radius is 2r + eps(), which guarantees that it is large enough.
    f32 r = 2 * bounding_radius + min_dist_to_edge + min_dist_to_vertex;

    v2* a = _DelaunayMeshAddVertex(mesh, r *cos(90 * DEG2RAD), r * sin(90 * DEG2RAD));
    v2* b = _DelaunayMeshAddVertex(mesh, r *cos(210 * DEG2RAD), r * sin(210 * DEG2RAD));
    v2* c = _DelaunayMeshAddVertex(mesh, r *cos(-30 * DEG2RAD), r * sin(-30 * DEG2RAD));

    struct QuarterEdge* ab = _DelaunayMeshAddEdge(mesh, a, b);
    struct QuarterEdge* bc = _DelaunayMeshAddEdge(mesh, b, c);
    struct QuarterEdge* ca = _DelaunayMeshAddEdge(mesh, c, a);

    QESplice(QESym(ab), bc);
    QESplice(QESym(bc), ca);
    QESplice(QESym(ca), ab);

    return mesh;
}

// ------------------------------------------------------------------------------------------------
void DeconstructDelaunayMesh(struct DelaunayMesh* mesh) {
    free(mesh->vertices);
    free(mesh->quarter_edges);
    free(mesh);
}

// ------------------------------------------------------------------------------------------------
size_t DelaunayMeshNumVertices(struct DelaunayMesh* mesh) { return mesh->n_vertices; }
size_t DelaunayMeshNumQuarterEdges(struct DelaunayMesh* mesh) { return mesh->n_quarter_edges; }
size_t DelaunayMeshNumEdges(struct DelaunayMesh* mesh) { return mesh->n_quarter_edges / 4; }

// ------------------------------------------------------------------------------------------------
const v2* DelaunayMeshGetVertex(struct DelaunayMesh* mesh, int i) { return mesh->vertices + i; }
struct QuarterEdge* DelaunayMeshGetQuarterEdge(struct DelaunayMesh* mesh, int i) { return mesh->quarter_edges + i; }

// ------------------------------------------------------------------------------------------------
struct QuarterEdge* DelaunayMeshGetQuarterEdgeBetween(struct DelaunayMesh* mesh, int i, int j) {
    if (i == j) {
        return NULL;  // Self edges do not exist
    }

    int n = DelaunayMeshNumVertices(mesh);
    if (i < 0 || i >= n || j < 0 || j >= n) {
        return NULL;  // Invalid index
    }

    // If we already have an edge between these two vertices, we are done.
    const v2* a = DelaunayMeshGetVertex(mesh, i);
    const v2* b = DelaunayMeshGetVertex(mesh, j);
    for (i = 0; i < mesh->n_quarter_edges; i++) {
        struct QuarterEdge* qe = mesh->quarter_edges + i;
        if (qe->vertex == a) {
            if (QESym(qe)->vertex == b) {
                return qe;
            }
        }
    }

    // Not found
    return NULL;
}

// ------------------------------------------------------------------------------------------------
bool DelaunayMeshHasEdge(struct DelaunayMesh* mesh, int i, int j) {
    return DelaunayMeshGetQuarterEdgeBetween(mesh, i, j) != NULL;
}

// ------------------------------------------------------------------------------------------------
v2* _DelaunayMeshAddVertex(struct DelaunayMesh* mesh, f32 x, f32 y) {
    if (mesh->n_vertices >= mesh->max_n_vertices) {
        return NULL; // At capacity
    }

    v2* vertex = &(mesh->vertices[mesh->n_vertices]);
    mesh->n_vertices++;

    vertex->x = x;
    vertex->y = y;
    return vertex;
}

// ------------------------------------------------------------------------------------------------
struct QuarterEdge* _DelaunayMeshAddEdge(struct DelaunayMesh* mesh, v2* a, v2* b) {
    if (mesh->n_quarter_edges >= mesh->max_n_quarter_edges - 3) {
        return NULL; // At capacity
    }

    struct QuarterEdge* ab = &(mesh->quarter_edges[mesh->n_quarter_edges]);
    struct QuarterEdge* lr = &(mesh->quarter_edges[mesh->n_quarter_edges + 1]);
    struct QuarterEdge* ba = &(mesh->quarter_edges[mesh->n_quarter_edges + 2]);
    struct QuarterEdge* rl = &(mesh->quarter_edges[mesh->n_quarter_edges + 3]);
    mesh->n_quarter_edges += 4;


    ab->vertex = a;
    ab->next = ab;
    ab->rot = lr;

    lr->next = rl;
    lr->rot = ba;
    
    ba->vertex = b;
    ba->next = ba;
    ba->rot = rl;
    
    rl->next = lr;
    rl->rot = ab;

    return ab;
}