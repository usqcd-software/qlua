#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "qlayout.h"                                                 /* DEPS */


static int eo_numsites(QDP_Lattice *lat, int node);

/* All lattice params are stored in mLattice.
 * Here we only keep a pointer to it.
 */
typedef struct {
    mLattice *S;
} params;

static void
get_lex_x(int *x, int l, int *s, int ndim)
{
    int i;

    for (i = ndim-1; i >= 0; --i) {
        x[i] = l % s[i];
        l = l / s[i];
    }
}

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
        for (i = 0; i < S->rank; i++) {
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

    int mc[QLUA_MAX_LATTICE_RANK];
    int i;

    S->node = QDP_this_node;
    node2coord(mc, QDP_this_node, S);

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

static int
eo_node_number(QDP_Lattice *lat, const int x[])
{
    params *p = QDP_get_lattice_params(lat);
    mLattice *S = p->S;
    int m[QLUA_MAX_LATTICE_RANK];
    int i, node;

    for (i = 0; i < S->rank; i++) {
        m[i] = ((x[i] + 1) * S->net[i] - 1) / S->dim[i];
    }
    node = coord2node(m, S);

    return node;
}

static int
eo_index(QDP_Lattice *lat, const int x[])
{
    params *p = QDP_get_lattice_params(lat);
    mLattice *S = p->S;
    int s = 0, l = 0, n = 1;
    int i;

    for (i = 0; i < S->rank; i++) {
        int m = ((x[i] + 1) * S->net[i] - 1) / S->dim[i];
        int x0 = (m * S->dim[i]) / S->net[i];
        int x1 = ((m + 1) * S->dim[i]) / S->net[i];
        l = l * (x1 - x0) + x[i] - x0;
        n = n * (x1 - x0);
        s += x[i];
    }

    if( s % 2==0 ) { /* even site */
        l /= 2;
    } else {
        l = (l + n) / 2;
    }
    return l;
}

static void
eo_get_coords(QDP_Lattice *lat, int x[], int node, int index)
{
    params *p = (params *) QDP_get_lattice_params(lat);
    mLattice *S = p->S;
    int nd = S->rank;
    int m[QLUA_MAX_LATTICE_RANK];
    int dx[QLUA_MAX_LATTICE_RANK];
    int sx[QLUA_MAX_LATTICE_RANK];
    int s0 = 0;
    int n0 = 1;
    int i;

    node2coord(m, node, S);
    
    for (i = 0; i < nd; i++) {
        x[i] = (m[i] * S->dim[i]) / S->net[i];
        int x1 = ((m[i] + 1) * S->dim[i]) / S->net[i];
        dx[i] = x1 - x[i];
        s0 += x[i];
        n0 *= dx[i];
    }

    int neven = (n0 + 1 - (s0 & 1)) / 2;
    if (index < neven) {
        int l = 2*index;
        int s1 = s0;

        get_lex_x(sx, l, dx, nd);
        for(i = 0; i<nd; ++i) s1 += sx[i];
        if ((s1&1)!=0) {
            get_lex_x(sx, l+1, dx, nd);
        }
    } else {
        int l = 2 * index - n0 + ((n0 & 1) * (s0 & 1));
        int s1 = s0;

        get_lex_x(sx, l, dx, nd);
        for (i = 0; i < nd; i++)
            s1 += sx[i];
        if ((s1&1)==0) {
            get_lex_x(sx, l+1, dx, nd);
        }
    }
    for (i = 0; i < nd; ++i) {
        x[i] += sx[i];
    }
}

void
qlua_sublattice(int lo[], int hi[], int node, void *env)
{
    mLattice *S = env;
    int nD[QLUA_MAX_LATTICE_RANK];
    int i;

    node2coord(nD, node, S);
    for (i = 0; i < S->rank; i++) {
        lo[i] = (nD[i] * S->dim[i]) / S->net[i];
        hi[i] = ((nD[i] + 1) * S->dim[i])/S->net[i];
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
