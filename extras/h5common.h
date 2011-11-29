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


/* attribute utilities */
int 
h5_dspace_equal(lua_State *L, hid_t s1, hid_t s2);
int     
h5_check_attr(lua_State *L, hid_t *p_a_id, hid_t obj_id, 
        const char *a_name, hid_t a_type, hid_t a_space);
int 
h5_check_attr_data(lua_State *L, hid_t *p_a_id, hid_t obj_id, 
        const char *a_name, const hid_t a_type, const hid_t a_space,
        const void *buf, hid_t mem_type);
char *
h5_make_str_data(lua_State *L, int n_str, int str_max_len, const char *s[]);
int
h5_check_attr_str_list(lua_State *L, hid_t *p_a_id, hid_t obj_id, 
        const char *a_name, int n_str, int str_max_len, const char *buf[]);
int
h5_check_attr_int(lua_State *L, hid_t *p_a_id, hid_t obj_id,
        const char *a_name, const int *int_buf);
int
h5_check_attr_array_int(lua_State *L, hid_t *p_a_id, hid_t obj_id, 
        const char *a_name, int ndims, const hsize_t dims[], const int *int_buf);
int
h5_check_attr_array1d_int(lua_State *L, hid_t *p_a_id, hid_t obj_id, 
        const char *a_name, int n1, const int *int_buf);
int
h5_check_attr_array2d_int(lua_State *L, hid_t *p_a_id, hid_t obj_id, 
        const char *a_name, int n1, int n2, const int *int_buf);


#endif/*H5COMMON_H_PLJNSRPMGIPXDQUJXFJE*/
