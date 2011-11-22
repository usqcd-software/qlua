#ifndef LAPH_COMMON_H_Trisdx0uvm7Huybtrmme
#define LAPH_COMMON_H_Trisdx0uvm7Huybtrmme

#include "qlua.h"                                                   /* DEPS */
#include "lattice.h"                                                /* DEPS */


#define LDIM 4      /* number of lattice dimensions */
#define NSPIN   4
#define NCOLOR  3


int 
calc_subgrid(QDP_Lattice *lat, int dims[], int cmin[]);

int 
qlua_check_laph_sol_list(lua_State *L, int tab_idx, mLattice *S,
               int *n_sol, int **tsrc, int **jvec, int **jspin,
               QDP_D3_DiracFermion ***sol, mLattice **have_S);

#endif/*LAPH_COMMON_H_Trisdx0uvm7Huybtrmme*/
