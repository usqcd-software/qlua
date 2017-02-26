#ifndef EXTRA_COMMON_H_ogugYvsijrsafuAfuaaw
#define EXTRA_COMMON_H_ogugYvsijrsafuAfuaaw

#include "modules.h"                                                /* DEPS */
#include "qlua.h"                                                   /* DEPS */
#include "lattice.h"                                                /* DEPS */

#define MAX_LDIM    8
int is_masternode();

int calc_subgrid(QDP_Lattice *lat, int dims[], int cmin[]);
int *locindex_space(QDP_Lattice *lat, int node, int t_axis, 
        const int locsize[], const int loc_c0[]);


void calc_D_C_eq_planewave(QDP_D_Complex *exp_ipx, const int mom[], const int x0[]);
void calc_F_C_eq_planewave(QDP_F_Complex *exp_ipx, const int mom[], const int x0[]);
#if ( QDP_Precision == 'F' )
# define calc_C_eq_planewave(a,b,c) calc_F_C_eq_planewave(a,b,c)
#elif ( QDP_Precision == 'D' )
# define calc_C_eq_planewave(a,b,c) calc_D_C_eq_planewave(a,b,c)
#else
# error Invalid QDP_Precision
#endif


#endif/*EXTRA_COMMON_H_ogugYvsijrsafuAfuaaw*/
