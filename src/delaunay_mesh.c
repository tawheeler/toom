#include "delaunay_mesh.h"

#include "constants.h"

#include <math.h>
#include <stdlib.h>

// ------------------------------------------------------------------------------------------------
f32 Det2x2(f32 a, f32 b, f32 c, f32 d) {
    return a * d - b * c;
}

// ------------------------------------------------------------------------------------------------
f32 Det3x3(f32 a, f32 b, f32 c,
        f32 d, f32 e, f32 f,
        f32 g, f32 h, f32 i) {
    return a * Det2x2(e, f, h, i) - b * Det2x2(d, f, g, i) + c * Det2x2(d, e, g, h);
}

// ------------------------------------------------------------------------------------------------
f32 Det4x4(f32 a, f32 b, f32 c, f32 d, 
           f32 e, f32 f, f32 g, f32 h,
           f32 i, f32 j, f32 k, f32 l,
           f32 m, f32 n, f32 o, f32 p) {
    return a * Det3x3(f, g, h,
                      j, k, l,
                      n, o, p)
         - b * Det3x3(e, g, h,
                      i, k, l,
                      m, o, p)
         + c * Det3x3(e, f, h,
                      i, j, l,
                      m, n, p)
         - d * Det3x3(e, f, g,
                      i, j, k,
                      m, n, o);
}

// ------------------------------------------------------------------------------------------------
f32 GetRightHandedness(const v2* a, const v2* b, const v2* c) {
    // return Det3x3(a->x, a->y, 1.0f, b->x, b->y, 1.0f, c->x, c->y, 1.0f);
    // return a->x * Det(b->y, 1.0f, c->y, 1.0f) - a->y * Det(b->x, 1.0f, c->x, 1.0f) + 1.0f * Det(b->x, b->y, c->x, c->y);
    // return a->x * Det(b->y, 1.0f, c->y, 1.0f) - a->y * Det(b->x, 1.0f, c->x, 1.0f) + Det(b->x, b->y, c->x, c->y);
    // return a->x * (b->y * 1.0f - 1.0f * c->y) - a->y * (b->x * 1.0f - 1.0f * c->x)  + (b->x * c->y - b->y * c->x)
    return a->x*(b->y - c->y) - a->y*(b->x - c->x) + (b->x*c->y - b->y*c->x);
}

// ------------------------------------------------------------------------------------------------
// Get the distance of point P to the line that goes through A and B.
// Note that A and B cannot be colocated.
f32 GetDistanceToLine(const v2* p, const v2* a, const v2* b) {
    v2 delta = sub(*b, *a);
    return abs(delta.x * (a->y - p->y) - delta.y * (a->x - p->x)) / length(delta);
}

// ------------------------------------------------------------------------------------------------
// Determine whether the given point is inside (or on the border of) the triangle defined by the
// three vertices.
// Returns 1 if inside the triangle, 0 if on its perimeter, and -1 if outside
int GetTriangleContainment(const v2* p, const v2* a, const v2* b, const v2* c) {
    f32 ab = GetRightHandedness(a, b, p);
    f32 bc = GetRightHandedness(b, c, p);
    f32 ca = GetRightHandedness(c, a, p);
    if (ab > 0 && bc > 0 && ca > 0) {
        return 1;
    } else if (ab >= 0 && bc >= 0 && ca >= 0) {
        return 0;
    } else {
        return -1;
    }
}

// ------------------------------------------------------------------------------------------------
// Determine whether the given point is inside the circle inscribed on the three vertices.
// Returns 1 if inside the circle,
//         0 if on the circle,
//    and -1 if outside the circle.
int GetCircleContainment(const v2* p, const v2* a, const v2* b, const v2* c) {
    float d = Det4x4(a->x, a->y, a->x*a->x + a->y*a->y, 1.0f,
                     b->x, b->y, b->x*b->x + b->y*b->y, 1.0f,
                     c->x, c->y, c->x*c->x + c->y*c->y, 1.0f,
                     p->x, p->y, p->x*p->x + p->y*p->y, 1.0f);

    if (d > 0) {
        return 1;
    } else if (d < 0) {
        return -1;
    } else {
        return 0;
    }
}

// ------------------------------------------------------------------------------------------------
bool IsDualEdge(const QuarterEdge* qe) { return qe->vertex == NULL; }

// ------------------------------------------------------------------------------------------------
bool IsPrimalEdge(const QuarterEdge* qe) { return qe->vertex != NULL; }

// ------------------------------------------------------------------------------------------------
QuarterEdge* QENext(const QuarterEdge* qe) { return qe->next; }
QuarterEdge* QERot(const QuarterEdge* qe) { return qe->rot; }
QuarterEdge* QESym(const QuarterEdge* qe) { return qe->rot->rot; }
QuarterEdge* QETor(const QuarterEdge* qe) { return qe->rot->rot->rot; }
QuarterEdge* QEPrev(const QuarterEdge* qe) { return qe->rot->next->rot; }
QuarterEdge* QELnext(const QuarterEdge* qe) { return QETor(qe)->next->rot; }

// ------------------------------------------------------------------------------------------------
void QESwapNexts(QuarterEdge* a, QuarterEdge* b) {
    QuarterEdge* tmp = a->next;
    a->next = b->next;
    b->next = tmp;
}

// ------------------------------------------------------------------------------------------------
void QESplice(QuarterEdge* a, QuarterEdge* b) {
    QESwapNexts(a->next->rot, b->next->rot);
    QESwapNexts(a, b);
}


// ------------------------------------------------------------------------------------------------
void QEFlipEdge(QuarterEdge* qe) {
    QuarterEdge* qe_sym = QESym(qe);
    QuarterEdge* qe_a = QEPrev(qe);
    QuarterEdge* qe_b = QEPrev(qe_sym);

    QESplice(qe, qe_a);
    QESplice(qe_sym, qe_b);
    QESplice(qe, QELnext(qe_a));
    QESplice(qe_sym, QELnext(qe_b));

    qe->vertex = QESym(qe_a)->vertex;
    qe_sym->vertex = QESym(qe_b)->vertex;
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
QuarterEdge* _DelaunayMeshAddEdge(struct DelaunayMesh* mesh, v2* a, v2* b) {
    if (mesh->n_quarter_edges >= mesh->max_n_quarter_edges - 3) {
        return NULL; // At capacity
    }

    QuarterEdge* ab = &(mesh->quarter_edges[mesh->n_quarter_edges]);
    QuarterEdge* lr = &(mesh->quarter_edges[mesh->n_quarter_edges + 1]);
    QuarterEdge* ba = &(mesh->quarter_edges[mesh->n_quarter_edges + 2]);
    QuarterEdge* rl = &(mesh->quarter_edges[mesh->n_quarter_edges + 3]);

    ab->index = mesh->n_quarter_edges;
    ab->vertex = a;
    ab->next = ab;
    ab->rot = lr;

    lr->index = mesh->n_quarter_edges + 1;
    lr->vertex = NULL;
    lr->next = rl;
    lr->rot = ba;
    
    ba->index = mesh->n_quarter_edges + 2;
    ba->vertex = b;
    ba->next = ba;
    ba->rot = rl;
    
    rl->index = mesh->n_quarter_edges + 3;
    rl->vertex = NULL;
    rl->next = lr;
    rl->rot = ab;

    mesh->n_quarter_edges += 4;

    return ab;
}

// ------------------------------------------------------------------------------------------------
struct DelaunayMesh* ConstructEmptyDelaunayMesh(
    f32 bounding_radius,
    f32 min_dist_to_vertex,
    f32 min_dist_to_edge,
    u32 max_n_vertices,
    u32 max_n_quarter_edges) {
    
    struct DelaunayMesh* mesh = (struct DelaunayMesh*) malloc(sizeof(struct DelaunayMesh));
    mesh->square_bounding_radius = bounding_radius*bounding_radius;
    mesh->square_min_dist_to_vertex = min_dist_to_vertex*min_dist_to_vertex;
    mesh->min_dist_to_edge = min_dist_to_edge;

    mesh->max_n_quarter_edges = max_n_quarter_edges;
    mesh->max_n_vertices = max_n_vertices;
    mesh->n_quarter_edges = 0;
    mesh->n_vertices = 0;

    mesh->vertices = (v2*) malloc(max_n_vertices * sizeof(v2));
    mesh->quarter_edges = (QuarterEdge*) malloc(max_n_quarter_edges * sizeof(QuarterEdge));

    // The triangle radius is 2r + eps(), which guarantees that it is large enough.
    f32 r = 2 * bounding_radius + min_dist_to_edge + min_dist_to_vertex;

    v2* a = _DelaunayMeshAddVertex(mesh, r *cos(90 * DEG2RAD), r * sin(90 * DEG2RAD));
    v2* b = _DelaunayMeshAddVertex(mesh, r *cos(210 * DEG2RAD), r * sin(210 * DEG2RAD));
    v2* c = _DelaunayMeshAddVertex(mesh, r *cos(-30 * DEG2RAD), r * sin(-30 * DEG2RAD));

    QuarterEdge* ab = _DelaunayMeshAddEdge(mesh, a, b);
    QuarterEdge* bc = _DelaunayMeshAddEdge(mesh, b, c);
    QuarterEdge* ca = _DelaunayMeshAddEdge(mesh, c, a);

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
QuarterEdge* DelaunayMeshGetQuarterEdge(struct DelaunayMesh* mesh, int i) { return mesh->quarter_edges + i; }

// ------------------------------------------------------------------------------------------------
QuarterEdge* DelaunayMeshGetQuarterEdgeBetween(struct DelaunayMesh* mesh, int i, int j) {
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
        QuarterEdge* qe = mesh->quarter_edges + i;
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
bool DelaunayMeshIsBoundaryVertex(const struct DelaunayMesh* mesh, const v2* v) {
    return v == mesh->vertices || v == mesh->vertices + 1 || v == mesh->vertices + 2;
}

// ------------------------------------------------------------------------------------------------
bool DelaunayMeshHasEdge(struct DelaunayMesh* mesh, int i, int j) {
    return DelaunayMeshGetQuarterEdgeBetween(mesh, i, j) != NULL;
}

// ------------------------------------------------------------------------------------------------
int DelaunayMeshAddVertex(struct DelaunayMesh* mesh, const v2* p) {
    // Ensure that the triangle is within bounds
    if (sqnorm(*p) >= mesh->square_bounding_radius) {
        return INVALID_MESH_INDEX;
    }

    // Get the enclosing triangle
    QuarterEdge* qe_dual = DelaunayMeshGetEnclosingTriangle2(mesh, p);
    QuarterEdge* qe_ab = qe_dual->rot;
    QuarterEdge* qe_bc = qe_dual->next->rot;
    QuarterEdge* qe_ca = qe_dual->next->next->rot;

    // Grab the vertices
    const v2 a = *(qe_ab->vertex);
    const v2 b = *(qe_bc->vertex);
    const v2 c = *(qe_ca->vertex);

    // If we are too close to an existing vertex, do nothing.
    if (min(sqnorm(sub(a, *p)), min(sqnorm(sub(b, *p)), sqnorm(sub(c, *p)))) < mesh->square_min_dist_to_vertex) {
        return INVALID_MESH_INDEX;
    }

    // Add the vertex
    v2* p_ptr = _DelaunayMeshAddVertex(mesh, p->x, p->y);

    // Check whether we are (effectively) on an edge.
    // If we are, we delete the existing edge and add four instead.
    f32 dist_to_ab = GetDistanceToLine(p, &a, &b);
    f32 dist_to_bc = GetDistanceToLine(p, &b, &c);
    f32 dist_to_ca = GetDistanceToLine(p, &c, &a);
    QuarterEdge* qe_start = NULL;
    if (min(dist_to_ab, min(dist_to_bc, dist_to_ca)) > mesh->min_dist_to_edge) {
        // Normal case. We are not too close to an edge.

        // Add the new edges and splice them in.
        QuarterEdge* ap = _DelaunayMeshAddEdge(mesh, qe_ab->vertex, p_ptr);
        QuarterEdge* bp = _DelaunayMeshAddEdge(mesh, qe_bc->vertex, p_ptr);
        QuarterEdge* cp = _DelaunayMeshAddEdge(mesh, qe_ca->vertex, p_ptr);

        QESplice(ap, qe_ab);
        QESplice(bp, qe_bc);
        QESplice(cp, qe_ca);

        // TODO: Figure out splice such that we get this desired effect.
        //       i.e. can we replace these next calls with three calls to splice?
        QuarterEdge* pa = QESym(ap);
        QuarterEdge* pb = QESym(bp);
        QuarterEdge* pc = QESym(cp);

        pa->next = bp->rot->rot;  // sym(bp)
        pb->next = cp->rot->rot;  // sym(bc)
        pc->next = ap->rot->rot;  // sym(ap)

        pa->rot->next = cp->rot;
        pb->rot->next = ap->rot;
        pc->rot->next = bp->rot;

        // Set our start quarter edge
        qe_start = pa;
    } else {
        // We are effectively on an edge.
        // Identify that edge, and cut it with the new vertex.

        // Let DE be the edge we are on, F be the far vertex qe, and G be the vertex qe across DE
        // from F.
        QuarterEdge* qe_dp = NULL;
        QuarterEdge* qe_ef = NULL;
        QuarterEdge* qe_fd = NULL;
        if (dist_to_ab < mesh->min_dist_to_edge) {
            qe_dp = qe_ab;
            qe_ef = qe_bc;
            qe_fd = qe_ca;
        } else if (dist_to_bc < mesh->min_dist_to_edge) {
            qe_dp = qe_bc;
            qe_ef = qe_ca;
            qe_fd = qe_ab;
        } else {
            qe_dp = qe_ca;
            qe_ef = qe_ab;
            qe_fd = qe_bc;
        }
        QuarterEdge* qe_ge = QEPrev(QESym(QEPrev(qe_dp)));

        // Our edge may not be the boundary edge
        int num_boundary_vertices =
            DelaunayMeshIsBoundaryVertex(mesh, qe_dp->vertex) + 
            DelaunayMeshIsBoundaryVertex(mesh, QESym(qe_dp)->vertex);
        if (num_boundary_vertices == 2) {
            return INVALID_MESH_INDEX;
        }

        // Grab another quarter edge we will need
        QuarterEdge* qe_dg = QEPrev(qe_dp);

        // Reset the things that point to DE
        qe_dg->next = QERot(qe_fd)->rot;  // sym(fd)
        QuarterEdge* qe_ge_tor = QETor(qe_ge);
        qe_ge_tor->next = QESym(qe_ef)->rot;
        qe_ef->next = QERot(qe_ge)->rot;  // sym(ge)
        QuarterEdge* qe_fd_tor = QETor(qe_fd);
        qe_fd_tor->next = QESym(qe_dg)->rot;

        // Reset PD as a quarter edge, and change it to DP
        QuarterEdge* qe_pd = QESym(qe_dp);
        qe_pd->vertex = p_ptr;
        {
            qe_pd->next = qe_pd;
            qe_pd->rot->next = qe_dp->rot;
            qe_dp->next = qe_dp;
            qe_dp->rot->next = qe_pd->rot;
        }

        // Create three new edges EP, FP, and F->P
        QuarterEdge* qe_ep = _DelaunayMeshAddEdge(mesh, qe_ef->vertex, p_ptr);
        QuarterEdge* qe_fp = _DelaunayMeshAddEdge(mesh, qe_fd->vertex, p_ptr);
        QuarterEdge* qe_gp = _DelaunayMeshAddEdge(mesh, qe_ge->vertex, p_ptr);

        // Splice them all correctly
        QESplice(qe_dp, qe_dg);
        QESplice(qe_gp, qe_ge);
        QESplice(qe_ep, qe_ef);
        QESplice(qe_fp, qe_fd);

        // TODO: Figure out splice such that we get this desired effect.
        //       i.e. can we replace these next calls with three calls to splice?
        qe_pd = QESym(qe_dp);
        QuarterEdge* qe_pe = QESym(qe_ep);
        QuarterEdge* qe_pf = QESym(qe_fp);
        QuarterEdge* qe_pg = QESym(qe_gp);

        qe_pd->next = QERot(qe_gp)->rot;  // sym(gp)
        qe_pe->next = QERot(qe_fp)->rot;  // sym(fp)
        qe_pf->next = QERot(qe_dp)->rot;  // sym(dp)
        qe_pg->next = QERot(qe_ep)->rot;  // sym(ep)

        QERot(qe_pd)->next = qe_fp->rot;
        QERot(qe_pe)->next = qe_gp->rot;
        QERot(qe_pf)->next = qe_ep->rot;
        QERot(qe_pg)->next = qe_dp->rot;

        // Set our start quarter edge
        qe_start = qe_pd;
    }

    // Check if we are locally Delaunay.
    // Walk around the outer edges and flip any edges that are not locally delaunay.
    // We can check an edge by seeing if the opposite vertex is within the inscribed circle
    // of the inner vertex + edge vertices.

    // We start at pa and rotate around it to get all edges that we have to check.
    // We need to walk the outer edges until we get back to pa.
    // Each outer edge is given by next(mesh, sym(mesh, qe))

    QuarterEdge* qe = qe_start;
    bool done = 0;
    while (!done) {
        QuarterEdge* qe_outer_edge = QESym(qe)->next;

        // Advance
        qe = qe->next;
        done = qe == qe_start;

        // Only consider the edge if it is not a bounding edge
        v2* src = qe_outer_edge->vertex;
        v2* dst = QESym(qe_outer_edge)->vertex;
        int num_boundary_vertices = 
            DelaunayMeshIsBoundaryVertex(mesh, src) + 
            DelaunayMeshIsBoundaryVertex(mesh, dst);
        if (num_boundary_vertices == 2) {  // one vertex being on the boundary is okay
            continue;
        }

        // Check the edge from qe_outer_edge to sym(qe_outer_edge)

        // Get the far vertex across the dividing edge
        v2* far = QESym(qe_outer_edge->next)->vertex;

        // If the edge contains a boundary vertex, don't flip it if it would
        // produce an inside-out triangle.
        if (num_boundary_vertices == 1) {
            if (GetTriangleContainment(dst, far, p, src) >= 0 ||
                GetTriangleContainment(dst, far, src, p) >= 0 ||
                GetTriangleContainment(src, far, dst, p) >= 0 ||
                GetTriangleContainment(src, far, p, dst) >= 0) {
                continue;
            }
        }

        if (GetCircleContainment(p, src, dst, far) > 0 ||
            GetCircleContainment(far, p, dst, src) > 0) {
            // Either p is inside the circle passing through src, dst, and far, or
            //  far is inside the circle passing through p, src, and dst.
            // We have to flip the edge.
            QEFlipEdge(qe_outer_edge);

            // We flipped the edge, so qe_outer_edge has to be traversed again.
            // Back it up.
            qe = QEPrev(qe);
            if (qe != qe_start) {
                qe = QEPrev(qe);
            }

            done = 0;
        }
    }

    // Return the index of the newly added vertex
    return mesh->n_vertices - 1;
}

// ------------------------------------------------------------------------------------------------
bool DelaunayMeshConstrainEdge(struct DelaunayMesh* mesh, int i, int j) {
    // TODO: If we have vertices intersecting the edge, we should split our constrained edge
    //       by those vertices and call this for each sub-edge.

    // TODO: If we try to constrain an edge that overlaps with an already-constrained edge,
    //       then we need to produce an intersection and fix it that way.

    if (i == j) {
        return 0;  // Cannot add self-edges
    }

    size_t n = DelaunayMeshNumVertices(mesh);
    if (i < 0 || i >= n || j < 0 || j >= n) {
        return 0;  // Invalid index
    }

    // If we already have an edge between these two vertices, we are done.
    v2* a = mesh->vertices + i;
    v2* b = mesh->vertices + j;
    QuarterEdge* qe_a = NULL;
    for (int qe_index = 0; qe_index < mesh->n_quarter_edges; qe_index++) {
        QuarterEdge* qe = mesh->quarter_edges + qe_index;
        if (qe->vertex == a) {
            qe_a = qe;
            if (QESym(qe)->vertex == b) {
                return 1;  // already exists
            }
        }
    }

    // --------------------------------------------------------------------------
    // The edge does not yet exist.
    // Walk around, from A, toward B and flip any offending edges.
    // Repeat until we have flipped our way to producing AB.

    for (int iter = 0; iter < 100; iter++) {
        
        // Rotate qe_a to the last coincident quarter edge that is CCW of the new segment.
        while (GetRightHandedness(a, QESym(qe_a)->vertex, b) <= 0.0) {
            qe_a = qe_a->next;
        }
        while (GetRightHandedness(a, QESym(qe_a->next)->vertex, b) > 0.0) {
            qe_a = qe_a->next;
        }

        // qe_a will then have a CCW triangle ACD.
        // If the far side of the triangle is B (ACD == ADB), then we are done.
        // Otherwise, D is on the other side of AB, so CD intersects AB, and we need to see if we
        // can flip CD.
        QuarterEdge* qe_dual =
            QETor(qe_a);  // The dual quarter edge that starts inside ACD and points across AC.
        v2* c = (qe_dual->next->rot)->vertex;
        v2* d = (qe_dual->next->next->rot)->vertex;

        if (d == b) {
            // We have B, so we have produced edge ab and are done.
            return 1;
        } else {
            // See if we can flip CD. We have to ensure that it is a convex quadrilateral.
            // E is the vertex on the far side.
            v2* e = QETor(QESym(qe_dual->next)->next)->vertex;
            if (GetRightHandedness(a, c, e) > 0 && GetRightHandedness(a, e, d) > 0) {
                // Flip it
                QuarterEdge* qe_cd = qe_dual->next->rot;
                QEFlipEdge(qe_cd);
            } else {
                // TODO: Handle this case. We need to progress past C and try the next triangle that
                //       could overlap.
                return 0;
            }
        }
    }

    // AAAAAAH this should never happen.
    return 0;
}

// ------------------------------------------------------------------------------------------------
QuarterEdge* DelaunayMeshGetEnclosingTriangle(
    struct DelaunayMesh* mesh,
    const v2* p,
    QuarterEdge* qe_dual) {

    // This should always terminate
    const int kMaxIters = 100;
    for (int iter = 0; iter < kMaxIters; iter++) {
        QuarterEdge* qe_ab = qe_dual->rot;
        QuarterEdge* qe_bc = qe_dual->next->rot;
        QuarterEdge* qe_ca = qe_dual->next->next->rot;

        const v2* a = qe_ab->vertex;
        const v2* b = qe_bc->vertex;
        const v2* c = qe_ca->vertex;

        if (GetRightHandedness(a, b, p) < 0) {
            qe_dual = QERot(qe_ab);  // Move across AB
        } else if (GetRightHandedness(b, c, p) < 0) {
            qe_dual = QERot(qe_bc);  // Move across BC
        } else if (GetRightHandedness(c, a, p) < 0) {
            qe_dual = QERot(qe_ca);  // Move across CA
        } else {
            return qe_dual;
        }
    }
    return qe_dual;
}

// ------------------------------------------------------------------------------------------------
QuarterEdge* DelaunayMeshGetEnclosingTriangle2(struct DelaunayMesh* mesh, const v2* p) {
    return DelaunayMeshGetEnclosingTriangle(mesh, p, QETor(mesh->quarter_edges));
}

// ------------------------------------------------------------------------------------------------
const v2* DelaunayMeshGetTriangleVertex1(struct DelaunayMesh* mesh, const QuarterEdge* qe_dual) {
    return (qe_dual->rot)->vertex;
}

// ------------------------------------------------------------------------------------------------
const v2* DelaunayMeshGetTriangleVertex2(struct DelaunayMesh* mesh, const QuarterEdge* qe_dual) {
    return (qe_dual->next->rot)->vertex;
}

// ------------------------------------------------------------------------------------------------
const v2* DelaunayMeshGetTriangleVertex3(struct DelaunayMesh* mesh, const QuarterEdge* qe_dual) {
    return (qe_dual->next->next->rot)->vertex;
}