
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

/* compare two dataspaces: return true (1) if they are equal
   both must be scalar or simple, have the same rank and dimensions
   otherwise, return false (0)
*/
int
h5_dspace_equal(lua_State *L, hid_t s1, hid_t s2) 
{
    if (!H5Sis_simple(s1) || !H5Sis_simple(s2))
        return 0;

    int r1 = H5Sget_simple_extent_ndims(s1),
        r2 = H5Sget_simple_extent_ndims(s2);
    if (r1 != r2)
        return 1;

    hsize_t *d1 = qlua_malloc(L, sizeof(hsize_t) * r1),
            *d2 = qlua_malloc(L, sizeof(hsize_t) * r2);
    H5Sget_simple_extent_dims(s1, d1, NULL);
    H5Sget_simple_extent_dims(s2, d2, NULL);

    int is_equal = 1;
    if (memcmp(d1, d2, sizeof(hsize_t) * r1)) 
        is_equal = 0;

    qlua_free(L, d1);
    qlua_free(L, d2);

    return is_equal;
}

/* Check that an attribute exists and check its datatype and dataspace,
   or, if it does not exist, create it;
   if p_a_id!=NULL: return attr handle, othewise close attr
    IN:
         obj_id     parent object (usually dataset)
         a_name     attr name
         a_type     datatype the attr must have
         a_space    dataspace the attr must have
        
    OUT:
         p_a_id     opened attr id 

    result:
       OK:
             1      attribute was created anew
             0      attribute existed and type/space are correct
       Fail: 
            -1      cannot open the attribute
            -2      wrong data type
            -3      wrong data space
*/
int
h5_check_attr(lua_State *L,
              hid_t *p_a_id, hid_t obj_id, const char *a_name, 
              hid_t a_type, hid_t a_space)
{
    hid_t a_id;
    int status = 0;

    if (H5Aexists(obj_id, a_name)) {
        a_id = H5Aopen(obj_id, a_name, H5P_DEFAULT);
        if (a_id < 0) {
            status = -1;
            goto clearerr_0;
        }
        if (! H5Tequal(a_type, H5Aget_type(a_id))) {
            status = -2;
            goto clearerr_1;
        }
        hid_t x_space   = H5Aget_space(a_id);
        if (! h5_dspace_equal(L, a_space, x_space)) {
            status = -3;
            goto clearerr_1;
        }
        status = 0;
    } else {
        a_id = H5Acreate(obj_id, a_name, a_type, a_space, H5P_DEFAULT, H5P_DEFAULT);
        if (a_id < 0) {
            status = -1;
            goto clearerr_0;
        }
        status = 1;
    }

    if (NULL != p_a_id)
        *p_a_id = a_id;
    else
        H5Aclose(a_id); /* XXX to avoid orphaned open attributes */

    return status;

clearerr_1:
    H5Aclose(a_id);
clearerr_0:
    if (NULL != p_a_id)
        *p_a_id = -1;

    return status;
}

/* Check that an attribute exists and check its datatype, dataspace, and data
   or, if it does not exist, create it and initialize with the data;
   if p_a_id!=NULL: return attr handle, othewise close attr

    IN:
         obj_id     parent object (usually dataset)
         a_name     attr name
         a_type     datatype the attr must have
         a_space    dataspace the attr must have
         buf        data the attr must have
         mem_type   data type of the attribute data in memory
        
    OUT:
         p_a_id     (ignored if ==NULL)
                    result=OK && p_a_id != NULL: opened attr id
                    otherwise: invalid hid_t; attribute is closed

    result:
       OK:
             1      attribute was created anew and initialized
             0      attribute existed and type, space, and data are correct
       Fail: 
            -1      cannot open+read or create+write the attribute
            -2      wrong data type
            -3      wrong data space
            -4      wrong data 
*/
int 
h5_check_attr_data(lua_State *L,
                   hid_t *p_a_id, hid_t obj_id, const char *a_name,
                   const hid_t a_type, const hid_t a_space,
                   void *buf, hid_t mem_type)
{
    hid_t a_id;
    int x_status = h5_check_attr(L, &a_id, obj_id, a_name, a_type, a_space);
    switch (x_status) {
        case -1:
        case -2:
        case -3:
            if (NULL != p_a_id)
                *p_a_id = -1;
            return x_status;
            break;
        case 0:
            {
                int buf_len = H5Tget_size(mem_type) * H5Sget_simple_extent_npoints(a_space);
                void *buf_r = qlua_malloc(L, buf_len);
                if (H5Aread(a_id, mem_type, buf_r)) {
                    x_status = -1;
                    qlua_free(L, buf_r);
                    goto clearerr_1;
                }
                if (memcmp(buf, buf_r, buf_len)) {
                    x_status = -4;
                    qlua_free(L, buf_r);
                    goto clearerr_1;
                }
                qlua_free(L, buf_r);
            } break;
        case 1:
            {
                if (H5Awrite(a_id, mem_type, buf)) {
                    x_status = -1;
                    goto clearerr_1;
                }
            } break;
        default: 
            /* FIXME unexpected result of check_attr; handle? */
            break;
    }

    if (NULL != p_a_id)
        *p_a_id = a_id;
    else
        H5Aclose(a_id); /* XXX to avoid orphaned open attributes */

    return x_status;

clearerr_1:
    H5Aclose(a_id);
    return x_status;
}

/* create string table : char [n_str][str_len_max] and populate it with 
   null-terminated strings from s
 */
char *
h5_make_str_data(lua_State *L, hsize_t ls_dims[2], const char *s[])
{
    int n_str = ls_dims[0], 
        str_len_max = ls_dims[1];
    char *d = qlua_malloc(L, sizeof(char) * n_str * str_len_max);
    assert(NULL != d);
    memset(d, '\0', n_str * str_len_max);
    for (int i = 0; i < n_str ; i++)
        strncpy(d + str_len_max * i, s[i], str_len_max);
    return d;
}
