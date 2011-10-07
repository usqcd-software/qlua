
#include "h5common.h"

#include <assert.h>
#include <string.h>

int 
is_hsizearray_equal(int n, hsize_t *a, hsize_t *b)
{
    for (; n-- ; a++, b++)
        if (*a != *b)
            return 0;
    return 1;
}


/* On the master node only,
   1) open a HDF5 file for writing or create if it does not exist
   2) create dset if it does not exist or check that the dimensions correspond 
      to required ones
   3) TODO add an option to use an opened file hid
 */
const char *
h5_open_write(lua_State *L,
        h5output *h5o,
        const char *h5file, const char *h5path,
        int rank, const hsize_t *dims)
{
    const char *err_str = NULL;

    H5E_auto_t old_ehandler;
    void *old_client_data;
//    hid_t error_stack = H5Eget_current_stack();
    hid_t error_stack = H5E_DEFAULT;

    if (0 != QDP_this_node)
        return NULL;

    assert(NULL != h5o &&
            NULL != h5file &&
            NULL != h5path);

    /* try opening the datafile */
    H5Eget_auto(error_stack, &old_ehandler, &old_client_data);
    H5Eset_auto(error_stack, NULL, NULL);
    h5o->h5f = H5Fopen(h5file, H5F_ACC_RDWR, H5P_DEFAULT);
    H5Eset_auto(error_stack, old_ehandler, old_client_data);

    if (h5o->h5f < 0) {    /* file does not exist */
        H5Eclear(error_stack);
        if (0 > (h5o->h5f = H5Fcreate(h5file, H5F_ACC_EXCL,
                            H5P_DEFAULT, H5P_DEFAULT))) 
            return "cannot create HDF5 file";
    } 

    h5o->rank       = rank;
    h5o->dims       = qlua_malloc(L, sizeof(h5o->dims[0]) * rank);
    memcpy(h5o->dims, dims, sizeof(h5o->dims[0]) * rank);

    if (0 > (h5o->dspace = H5Screate_simple(h5o->rank, h5o->dims, NULL))) {
        err_str = "cannot create dataspace" ;
        goto clearerr_1;
    }

    /* try opening the dataset */
    H5Eget_auto(error_stack, &old_ehandler, &old_client_data);
    H5Eset_auto(error_stack, NULL, NULL);
    h5o->dset = H5Dopen(h5o->h5f, h5path, H5P_DEFAULT);
    H5Eset_auto(error_stack, old_ehandler, old_client_data);

    if (0 < h5o->dset) { /* dataset exists */
        hid_t dspace_orig;
        if (0 > (dspace_orig = H5Dget_space(h5o->dset))) {
            err_str = "cannot get dataspace from dataset";
            goto clearerr_1;
        }
        int rank_orig;
        rank_orig = H5Sget_simple_extent_ndims(dspace_orig);
        if (rank != rank_orig) {
            err_str = "dspace rank is not consistent";
            goto clearerr_1;
        }
        hsize_t *dims_orig = NULL;
        dims_orig = qlua_malloc(L, sizeof(dims_orig[0]) * rank_orig);
        H5Sget_simple_extent_dims(dspace_orig, dims_orig, NULL);
        if (!is_hsizearray_equal(rank_orig, dims_orig, h5o->dims)) {
            err_str = "dspace dims are not consistent";
            qlua_free(L, dims_orig);
            goto clearerr_1;
        }
        qlua_free(L, dims_orig);
    } else { /* dataset does not exist */
        H5Eclear(error_stack);
        if (0 > (h5o->dset = H5Dcreate(h5o->h5f, h5path, H5T_IEEE_F64BE, 
                        h5o->dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT))) {
            err_str = "cannot create dataset" ;
            goto clearerr_1;
        }
    }

    return NULL;

clearerr_1:
    qlua_free(L, h5o->dims);
    return err_str;
}
/* close HDF5 structures 
   XXX h5o is NOT deallocated
 */
const char *
h5_close(lua_State *L, h5output *h5o)
{
    if (0 != QDP_this_node)
        return NULL;

    qlua_free(L, h5o->dims);
    h5o->dims = NULL;

    if (H5Sclose(h5o->dspace) 
            || H5Dclose(h5o->dset)
            || H5Fclose(h5o->h5f))
        return "cannot close HDF5 file, dataspace or dataset";
    h5o->dspace = h5o->dset = h5o->h5f = -1;
    return NULL;
}
