#include "modules.h"                                                /* DEPS */
#include "qlua.h"                                                   /* DEPS */
#include "lattice.h"                                                /* DEPS */
#include "extras_common.h"

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
