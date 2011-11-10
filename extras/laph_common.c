

#include "laph_common.h"
#include <assert.h>

/* calculate subgrid dimensions and "lowest corner" coordinate for `lat'
   TODO check that the subgrid is rectilinear 
 */
int 
calc_subgrid(QDP_Lattice *lat, int dims[], int cmin[])
{
    assert(LDIM == QDP_ndim_L(lat));
    int n_site = QDP_sites_on_node_L(lat);
    assert(0 < n_site);

    QDP_get_coords_L(lat, cmin, QDP_this_node, 0);
    int cmax[LDIM];
    for (int k = 0 ; k < LDIM ; k++)
        cmax[k] = cmin[k];
    
    int coord[LDIM];
    for (int i_site = 1; i_site < n_site ; i_site++) {
        QDP_get_coords_L(lat, coord, QDP_this_node, i_site);
        for (int k = 0 ; k < LDIM ; k++) {
            if (coord[k] < cmin[k])
                cmin[k] = coord[k];
            if (cmax[k] < coord[k])
                cmax[k] = coord[k];
        }
    }
    
    int vol = 1;
    for (int k = 0; k < LDIM; k++) {
        dims[k] = 1 + cmax[k] - cmin[k];
        assert(0 < dims[k]);
        vol *= dims[k];
    }
    assert(vol == n_site);

    return 0;
}
