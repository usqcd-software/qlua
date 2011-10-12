#include "qlua.h"                                                    /* DEPS */
#include "qdp.h"
#include "qmp.h"
#include "lattice.h"                                                 /* DEPS */
#include "qlayout.h"                                                 /* DEPS */

extern void XXXdump(const char *fmt, ...);

/* All lattice params are stored in mLattice.
 * Here we only keep a pointer to it.
 */
typedef struct {
    mLattice *S;
    int numsites;
} params;

#if 0 /* XXXX */
static void
get_lex_x(int *x, int l, int *s, int ndim)
{
    int i;
    for (i = ndim-1; i >= 0; --i) {
        x[i] = l % s[i];
        l = l / s[i];
    }
}
#endif /* XXX */

static void
node2coord(int *x, int n, mLattice *S)
{
    int i;

    for(i = 0; i < S->rank; i++) {
        x[i] = n % S->net[i];
        n = n / S->net[i];
    }
}

static int
coord2node(int *x, mLattice *S)
{
    int i;
    int l = 0;
    int f = 1;
    for (i = 0; i < S->rank; i++) {
        l = l + f * x[i];
        f = f * S->net[i];
    }
    return l;
}

static int prime[] = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71
};
#define MAXPRIMES (sizeof(prime)/sizeof(int))

static void
eo_setup(QDP_Lattice *lat, void *args)
{
    mLattice *S = args;
    QDP_allocate_lattice_params(lat, sizeof (params));
    params *p = QDP_get_lattice_params(lat);

    p->S = S;

    if (QMP_get_msg_passing_type() != QMP_SWITCH) {
        int nd2 = QMP_get_allocated_number_of_dimensions();
        const int *nsquares2 = QMP_get_allocated_dimensions();
        int i;

        for (i = 0; i < S->rank; i++) {
            S->net[i] = (i < nd2) ? nsquares2[i] : 1;
        }
    } else { /* not QMP_GRID */
        int squaresize[QLUA_MAX_LATTICE_RANK];
        int extrafactors[QLUA_MAX_LATTICE_RANK];
        int i;

        for (i=0; i < S->rank; i++) {
            squaresize[i] = S->dim[i];
            extrafactors[i] = 1;
            S->net[i] = 1;
        }

        /* Figure out dimensions of rectangle */
        int n = QMP_get_number_of_nodes();   /* nodes to factor */
        int k = MAXPRIMES-1;
        while (n > 1) {
            /* figure out which prime to divide by starting with largest */
            /* if no factor found, assume n is prime */
            while ((k >= 0) && (n % prime[k] != 0)) --k;
            int pfac = (k>=0) ? prime[k] : n;

            /* figure out which direction to divide */
            /* find largest divisible dimension of h-cubes */
            /* if one direction with largest dimension has already been
               divided, divide it again.  Otherwise divide first direction
               with largest dimension. */
            int j = -1;
            int i;

            for (i = 0; i < S->rank; i++) {
                if (squaresize[i] % pfac == 0) {
                    if ((j<0) ||
                        (extrafactors[j] * squaresize[i] > 
                         extrafactors[i] * squaresize[j])) {
                        j = i;
                    } else if (extrafactors[j] * squaresize[i] == 
                               extrafactors[i] * squaresize[j]) {
                        if ((S->net[j] == 1) || (S->net[i] != 1))
                            j = i;
                    }
                }
            }

            /* This can fail if we run out of prime factors in the dimensions */
            /* then just choose largest dimension */
            if (j < 0) {
                int i;

                for (i = 0; i < S->rank; i++) {
                    if ((j<0) ||
                        (extrafactors[j] * squaresize[i] >
                         extrafactors[i] * squaresize[j]) ) {
                        j = i;
                    } else if (extrafactors[j] * squaresize[i] ==
                               extrafactors[i] * squaresize[j]) {
                        if((S->net[j] == 1) || (S->net[i] != 1))
                            j = i;
                    }
                }
                n /= pfac;
                extrafactors[j] *= pfac;
                S->net[j] *= pfac;
            } else {
                n /= pfac;
                squaresize[j] /= pfac;
                S->net[j] *= pfac;
            }
        }
    } /* not QMP_GRID */

    int numsites = 1;
    int mc[QLUA_MAX_LATTICE_RANK];
    int i;

    node2coord(mc, QDP_this_node, S);
    for (i=0; i < S->rank; i++) {
        int x0 = (mc[i] * S->dim[i] + S->net[i] - 1) / S->net[i];
        int x1 = ((mc[i]+1) * S->dim[i] + S->net[i] - 1)/S->net[i];
        numsites *= x1 - x0;
    }
    p->numsites = numsites;
    S->node = QDP_this_node;
    for (i = 0; i < S->rank; i++) {
        int x = mc[i];

        mc[i] = x + 1;
        if (mc[i] == S->net[i])
            mc[i] = 0;
        S->neighbor_up[i] = coord2node(mc, S);

        mc[i] = x - 1;
        if (mc[i] < 0)
            mc[i] = S->net[i] - 1;
        S->neighbor_down[i] = coord2node(mc, S);

        mc[i] = x;
    }
}

static void
eo_free(QDP_Lattice *lat)
{
}

static int
eo_numsites(QDP_Lattice *lat, int node)
{
    params *p = QDP_get_lattice_params(lat);
    mLattice *S = p->S;

    if (node == QDP_this_node) {
        return p->numsites;
    } else {
        int numsites = 1;
        int nd = S->rank;
        int mc[QLUA_MAX_LATTICE_RANK];
        int i;

        node2coord(mc, node, S);
        for (i = 0; i<nd; ++i) {
            int x0 = (mc[i] * S->dim[i]) / S->net[i];
            int x1 = ((mc[i] + 1) * S->dim[i]) / S->net[i];
            numsites *= x1-x0;
        }
        return numsites;
    }
}

static int
eo_node_number(QDP_Lattice *lat, const int x[])
{
    params *p = QDP_get_lattice_params(lat);
    mLattice *S = p->S;
    int m[QLUA_MAX_LATTICE_RANK];
    int i, n;
    
    for (i = 0; i < S->rank; i++) {
        m[i] = (x[i] * S->net[i]) / S->dim[i];
    }
    n = coord2node(m, S);

    return n;
}

/* Used by clover and mdwf as well */
void
qlua_sublattice(int lo[], int hi[], int node, void *env)
{
    mLattice *S = env;
    int nD[QLUA_MAX_LATTICE_RANK];
    int i;

    node2coord(nD, node, S);
    for (i = 0; i < S->rank; i++) {
        lo[i] = (nD[i] * S->dim[i] + S->net[i] - 1) / S->net[i];
        hi[i] = ((nD[i] + 1) * S->dim[i] + S->net[i] - 1)/S->net[i];
    }
}

static int
eo_index(QDP_Lattice *lat, const int x[])
{
    params *lp = QDP_get_lattice_params(lat);
    mLattice *S = lp->S;
    int s = 0, l = 0, z  = 0, f = 1;
    int i, p;

    for (i = 0; i < S->rank; i++) {
        int n = (x[i] * S->net[i]) / S->dim[i];
        int lo = (n * S->dim[i] + S->net[i] - 1) / S->net[i];
        int hi = ((n + 1) * S->dim[i] + S->net[i] - 1)/S->net[i];
        z += lo;
        s += x[i];
        l = l + (x[i] - lo) * f;
        f *= hi - lo;
    }
    z &= 1;
    s &= 1;

    p = l / 2 + s * (f + 1 - z) / 2;

    return p;
}

static void
eo_get_coords(QDP_Lattice *lat, int x[], int node, int index)
{
    params *lp = QDP_get_lattice_params(lat);
    mLattice *S = lp->S;
    int lo[QLUA_MAX_LATTICE_RANK];
    int hi[QLUA_MAX_LATTICE_RANK];
    int s = 0, z = 0, f = 1;
    int i, e, l;

    qlua_sublattice(lo, hi, node, S);
    for (i = 0; i < S->rank; i++) {
        z += lo[i];
        f *= hi[i] - lo[i];
    }
    z &= 1;
    e = (f + 1 - z) / 2;
    if (index < e) {
        s = 0;
    } else {
        s = 1;
        index -= e;
    }
    l = s + z + 2 * (index - s * z);
    for (i = 0; i < S->rank; i++) {
        int d = hi[i] - lo[i];
        x[i] = lo[i] + l % d;
        l = l / d;
    }
}

QDP_Layout qlua_layout = {
    eo_setup,
    eo_free,
    eo_numsites,
    eo_node_number,
    eo_index,
    eo_get_coords
};
