#include "modules.h"                                                /* DEPS */
#include "qlua.h"                                                   /* DEPS */
#include "lattice.h"                                                /* DEPS */
#include "extras_common.h"

#include <assert.h>

int
is_masternode()
{
    return (0 == QDP_this_node ||       /* masternode */
            QDP_this_node < 0);         /* qdp not initialized */
}

/*****************************************************************************
 coordinate functions
 *****************************************************************************/

int 
calc_subgrid(QDP_Lattice *lat, int dims[], int cmin[])
/* calculate subgrid dimensions and "lowest corner" coordinate for `lat'
   TODO check that the subgrid is rectilinear 
 */
{
    int ndim = QDP_ndim_L(lat);
    assert(ndim <= MAX_LDIM);
    int n_site = QDP_sites_on_node_L(lat);
    assert(0 < n_site);

    QDP_get_coords_L(lat, cmin, QDP_this_node, 0);
    int cmax[MAX_LDIM];
    for (int k = 0 ; k < ndim ; k++)
        cmax[k] = cmin[k];
    
    int coord[MAX_LDIM];
    for (int i_site = 1; i_site < n_site ; i_site++) {
        QDP_get_coords_L(lat, coord, QDP_this_node, i_site);
        for (int k = 0 ; k < ndim ; k++) {
            if (coord[k] < cmin[k])
                cmin[k] = coord[k];
            if (cmax[k] < coord[k])
                cmax[k] = coord[k];
        }
    }
    
    int vol = 1;
    for (int k = 0; k < ndim; k++) {
        dims[k] = 1 + cmax[k] - cmin[k];
        assert(0 < dims[k]);
        vol *= dims[k];
    }
    assert(vol == n_site);

    return 0;
}

int *
locindex_space(QDP_Lattice *lat, int node, int t_axis, 
        const int locsize[], const int loc_c0[])
/* calculate local index -> local space index assuming rectangular subgrid
    lat         QDP lattice object
    ndim        == lat. dimension D (sizes of coord, c0, latsize)
    t_axis      skipped dimension; if 0 < t_axis or ndim<=t_axis, full index is calculated
    locsize     local lattice size
    loc_c0      coor. D-vector: min.corner of the local subgrid 
    locsize, loc_c0 are computed by calc_subgrid
    RETURN : pointer to the index or NULL if ENOMEM
 */
{
    int ndim = QDP_ndim_L(lat);
    assert(ndim <= MAX_LDIM);
    int dim_stride[MAX_LDIM],
        coord[MAX_LDIM];
    int fact;
    /* calculate strides for dimensions */
    for (int d = 0, fact = 1 ; d < ndim ; d++) {
        if (d == t_axis)
            dim_stride[d]   = 0;
        else {
            dim_stride[d]   = fact;
            fact            *= locsize[d];
        }
    }
    int n_site  = QDP_numsites_L(lat, QDP_this_node);
    int *idx    = malloc(n_site * sizeof(int));
    if (NULL == idx)
        return NULL;
    for (int i_site = 0 ; i_site < n_site ; i_site++) {
        QDP_get_coords_L(lat, coord, node, i_site);
        int i_x3 = 0;
        for (int d = 0 ; d < ndim ; d++)
            i_x3 += dim_stride[d] * (coord[d] - loc_c0[d]);
        idx[i_site] = i_x3;
    }
    return idx;
}



struct get_exp_ipx_C_arg {
    int dim;
    const int *latsize;
    const int *ft_mom;
    const int *ft_x0;
};

static void
get_planewave_C_D(QLA_D_Complex *dest, int coord[], void *arg_)
{
    struct get_exp_ipx_C_arg *arg = (struct get_exp_ipx_C_arg *)arg_;
    double phase = 0.;
    for (int d = 0; d < arg->dim ; d++) {
        phase += 2 * M_PI * (coord[d] - arg->ft_x0[d]) * 
                (double)arg->ft_mom[d] / arg->latsize[d];
    }
    QLA_c_eq_r_plus_ir(*dest, cos(phase), sin(phase));
}

void
get_planewave_C_F(QLA_F_Complex *dest, int coord[], void *arg_)
{
    struct get_exp_ipx_C_arg *arg = (struct get_exp_ipx_C_arg *)arg_;
    double phase = 0.;
    for (int d = 0; d < arg->dim ; d++) {
        phase += 2 * M_PI * (coord[d] - arg->ft_x0[d]) * 
                (double)arg->ft_mom[d] / arg->latsize[d];
    }
    QLA_c_eq_r_plus_ir(*dest, cos(phase), sin(phase));
}

void
calc_F_C_eq_planewave(QDP_F_Complex *exp_ipx, 
        const int mom[], const int x0[])
{
    QDP_Lattice *lat = QDP_F_get_lattice_C(exp_ipx);
    int dim = QDP_ndim_L(lat);
    int *latsize = malloc(sizeof(int) * dim);
    QDP_latsize_L(lat, latsize);
    struct get_exp_ipx_C_arg arg = { dim, latsize, mom, x0 };
    QDP_F_C_eq_funca(exp_ipx, get_planewave_C_F, &arg, QDP_all_L(lat));
    free(latsize);
}
void
calc_D_C_eq_planewave(QDP_D_Complex *exp_ipx, 
        const int mom[], const int x0[])
{
    QDP_Lattice *lat = QDP_D_get_lattice_C(exp_ipx);
    int dim = QDP_ndim_L(lat);
    int *latsize = malloc(sizeof(int) * dim);
    QDP_latsize_L(lat, latsize);
    struct get_exp_ipx_C_arg arg = { dim, latsize, mom, x0 };
    QDP_D_C_eq_funca(exp_ipx, get_planewave_C_D, &arg, QDP_all_L(lat));
    free(latsize);
}
