#ifndef H5COMMON_H_PLJNSRPMGIPXDQUJXFJE
#define H5COMMON_H_PLJNSRPMGIPXDQUJXFJE

#include <hdf5.h>
#include "qlua.h"

typedef struct {
    hid_t       h5f;
    hid_t       dset;
    
    int         rank;
    hsize_t     *dims;
    hid_t       dspace;

} h5output;

int       
is_hsizearray_equal(int n, hsize_t *a, hsize_t *b);

const char *
h5_open_write(lua_State *L, 
        h5output *h5o,      
        const char *h5file, const char *h5path,
        int rank, const hsize_t *dims);

const char *
h5_close(lua_State *L, h5output *h5o);

#endif/*H5COMMON_H_PLJNSRPMGIPXDQUJXFJE*/
