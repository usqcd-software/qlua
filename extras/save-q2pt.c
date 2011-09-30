
/* TODO separate h5output code from squark and save-q2pt */

typedef struct {
    hid_t       h5f;
    hid_t       dset;
    
    int         rank;
    hsize_t     *dims;
    hid_t       dspace;

} h5output;


/* On the master node only,
   1) open a HDF5 file for writing or create if it does not exist
   2) create dset if it does not exist or check that the dimensions correspond 
      to required ones
   3) TODO add an option to use an opened file hid
 */
static const char *
h5_open_write(luaL_state *L, 
        h5output *h5o,
        const char *h5file, const char *h5path,
        int rank, const hsize_t *dims)
{
    const char *err_str = NULL;

    H5E_auto_t old_ehandler;
    void *old_client_data;
    hid_t error_stack = H5Eget_current_stack();

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
        if (6 != rank_orig) {
            err_str = "dspace rank is not consistent";
            goto clearerr_1;
        }
        hsize_t dims_orig[6];
        H5Sget_simple_extent_dims(dspace_orig, dims_orig, NULL);
        if (!is_hsizearray_equal(rank_orig, dims_orig, h5o->dims)) {
            err_str = "dspace dims are not consistent";
            goto clearerr_1;
        }
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
static const char *
h5_close(h5output *h5o)
{
    if (0 != QDP_this_node)
        return NULL;

    qlua_free(L, h5o->dims);

    if (H5Sclose(h5o->dspace) 
            || H5Dclose(h5o->dset)
            || H5Fclose(h5o->h5f))
        return "cannot close HDF5 file, dataspace or dataset";
    else
        return NULL;
}



#define NS      4
#define LDIM    4

typedef struct {
    h5output    h5o;

    int lt;
    int n_v;
    int i_vec_prop;
    int i_spin_prop;

    int         buf_rank;
    hsize_t     buf_dims[4];
    hid_t       buf_dspace;
} q2pt_h5output;

static const char *
q2pt_h5_open(luaL_state *L,
            qlua_h5output *q2pt_h5o,
            const char *h5file, const char *h5path,
            int lt, int n_v, int i_vec_prop, int i_spin_prop)
{
    const char *err_str = NULL;

    hsize_t dims[7] = { lt, lt, n_v, n_v, NS, NS, 2};
    err_str = h5_open_write(L, &(q2pt_h5o->h5o), h5file, h5path, 7, dims);
    if (NULL != err_str)
        return err_str;

    q2pt_h5o->buf_rank = 4;
    q2pt_h5o->buf_dims[0]   = lt;
    q2pt_h5o->buf_dims[1]   = n_v;
    q2pt_h5o->buf_dims[2]   = NS;
    q2pt_h5o->buf_dims[3]   = 2;

    q2pt_h5o->i_vec_prop    = i_vec_prop;
    q2pt_h5o->i_spin_prop   = i_spin_prop;

    if (0 > (q2pt_h5o->buf_dspace = 
                H5Screate_simple(q2pt_h5o->buf_rank, q2pt_h5o->buf_dims, NULL)))
        return "cannot create memory dataspace";

    return NULL;
}


/* close q2pt_h5 structures
   XXX h5o is NOT deallocated
 */
static const char *
q2pt_h5_close(luaL_state *L, qlua_h5output *q2pt_h5o)
{
    if (H5Sclose(q2pt_h5o->buf_dspace))
        return "cannot close HDF5 membuf";
    return h5_close(L, &(q2pt->h5o));
}

/* write a portion of data */    
static const char *
q2pt_h5_write(luaL_state *L,
             qlua_h5output *q2pt_h5o, )
{
    /* TODO */
}


/* main:
   cycle through propagator(s) components and save its projection(s) on 
   LapH eigenmodes for each timeslice
 */
static const char *
save_q2pt(lua_State *L, 
        mLattice *S,
        const char *h5_file,
        const char *h5_path,
        int t0_prop,
        int i_v_prop, int i_s_prop, QDP_D3_DiracFermion *sol,
        int n_v, QDP_D3_ColorVector **v,
        int t_axis)
{
    /* TODO */
}

/* QLua wrapper */
int         
q_save_q2pt(lua_State *L)
{
    /* TODO */
}
