/* Only Nc=3 reader is supported */
#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "latcolmat.h"                                               /* DEPS */
#include "qend.h"                                                    /* DEPS */
#include "qmp.h"
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>

#if USE_Nc3
static const char milc_io[] = "milc";

static int
milcio_error(lua_State *L, const char *errmsg)
{
    return luaL_error(L, "L:read_MILC_gauge() error: %s", errmsg);
}


/* FIXME duplicated many times over(nersc_io.c, ???) */
static void
site2coord(int *coord, long long site, int nd, const int *dim)
{
    int i;

    for (i = 0; i < nd; i++) {
        coord[i] = site % dim[i];
        site = site / dim[i];
    }
}

static void
milc_cksum_update(uint32_t *cksum29, uint32_t *cksum31, uint32_t *data, size_t len, size_t offset)
{
    for (int i = 0 ; i < len ; i++) {
        size_t k = offset + i;
        int k29 = k % 29,
            k31 = k % 31;
        uint32_t x = data[i];
        *cksum29    ^= (x << k29) + (x >> (32 - k29));
        *cksum31    ^= (x << k31) + (x >> (32 - k31));
    }
}

/* 
   if broadcast str == NULL, then len=0 is communicated and NULL returned on slaves
   empty string always returns with len=1
 */
static char *
broadcast_string(lua_State *L, const char *str)
{
    int len = (NULL != str ? 1 + strlen(str) : 0);
    QMP_broadcast(&len, sizeof(len));
    if (0 == len)
        return NULL;

    char *s1 = NULL;
    if (NULL == (s1 = qlua_malloc(L, len))) {
        luaL_error(L, "broadcast_string: memory allocation error");
        return NULL;
    }
    if (NULL != str) {
        strncpy(s1, str, len);
        s1[len - 1] = '\0';
    }
    else s1[0] = '\0';
    QMP_broadcast(s1, len);
    return s1;
}

static int
read_milc_gauge_master(
                lua_State *L, 
                mLattice *S,
                char **datetmp1,
                QLA_D3_ColorMatrix **U,
                const char *name)
{
    char *status = NULL,
         *stat_b = NULL;
    char *datetmp = NULL;

    if (4 != S->rank) 
        return milcio_error(L, "not implemented for dim!=4");

    int milc_nc = 3;
    int  hdrbuf_len = 96;
    uint32_t *f_dim, 
             cksum29, f_cksum29,
             cksum31, f_cksum31;
    int byterev;

    int site_len = S->rank * 2*milc_nc*milc_nc,
        real_size = sizeof(float) ;
    char *site_buf = NULL;

    char *key, *value;
    int  *coord = NULL;

    long long volume;
    long long site;
    int s_node;
    int i;
    QLA_D3_ColorMatrix *CM = NULL;

    QLA_D_Complex cone;
    QLA_c_eq_r_plus_ir(cone, 1.0, 0.0);
    QLA_D3_ColorMatrix mone;
    QLA_D3_M_eq_c(&mone, &cone);

    /* open file and read header */
    FILE *f = fopen(name, "rb");
    if (f == 0) {
        status = "file open error";
        goto clearerr_0;
    }

    struct __attribute__((__packed__)) {
        uint32_t magic;
        uint32_t f_dim[4];
        char datetmp[64];
        uint32_t order;
        uint32_t cksum[2];
    } hdrbuf;

    if (1 != fread(&hdrbuf, sizeof(hdrbuf), 1, f)) {
        status = "missing header";
        goto clearerr_0;
    }
    byterev = 0;

    if (20103 != hdrbuf.magic) {
        swap_endian((char *)&hdrbuf.magic, 4, 1);  
        if (20103 != hdrbuf.magic) {
            status = "bad magic number in MILC gauge file";
            goto clearerr_0;
        }
        byterev = 1;
    }
    if (byterev) {
        swap_endian((char *)&hdrbuf.f_dim, 4, 4);
        swap_endian((char *)&hdrbuf.order, sizeof(hdrbuf.order), 1);
        swap_endian((char *)&hdrbuf.cksum, sizeof(hdrbuf.cksum[0]), 2);
    }       
    /* check dimensions */
    for (i = 0 ; i < S->rank ; i++) {
        if (S->dim[i] != hdrbuf.f_dim[i]) {
            status = "dimension mismatch";
            goto clearerr_0;
        }
    }
    /* copy date */
    datetmp = (char *) qlua_malloc(L, 65);
    memcpy(datetmp, hdrbuf.datetmp, 64);
    datetmp[64] = '\0';

    /* Read the data and send it to the target host */
    for (volume = 1, i = 0; i < S->rank; i++)
        volume *= S->dim[i] ;
    /* read every site in order on the master
     * compute the checksum and send it to the target node
     */
    coord = qlua_malloc(L, S->rank * sizeof (int));
    site_buf = qlua_malloc(L, site_len*real_size);
    CM = qlua_malloc(L, S->rank * sizeof (QLA_D3_ColorMatrix));
    QMP_msgmem_t mm = QMP_declare_msgmem(&CM[0], S->rank * sizeof (CM[0]));

    /* Go through all sites */
    cksum29 = 0;
    cksum31 = 0;
    for (site = 0; site < volume; site++) {
        if ((status == NULL) && (1 != fread(site_buf, site_len * real_size, 1, f))) {
            status = "file read error";
            goto clearerr_1;
        }
        if (byterev)
            swap_endian(site_buf, real_size, site_len);
        /* collect the checksum */
        milc_cksum_update(&cksum29, &cksum31, (uint32_t *)site_buf, site_len, site * site_len);

        for (int i = 0 ; i < S->rank ; i++) {
            for (int a = 0 ; a < milc_nc ; a++)
                for (int b = 0 ; b < milc_nc ; b++) {
                    /* FIXME is it C-matr? */
                    float *r1 = (float *)site_buf + 2*(b + milc_nc*(a + milc_nc* i));
                    QLA_c_eq_r_plus_ir(QLA_elem_M(CM[i], a, b), r1[0], r1[1]);
                }
        }

        /* place ColorMatrix in U where it belongs */
        site2coord(coord, site, S->rank, S->dim);
        s_node = QDP_node_number_L(S->lat, coord);
        if (s_node == QDP_this_node) {
            int idx = QDP_index_L(S->lat, coord);
            int d;
            
            for (d = 0; d < S->rank; d++)
                QLA_D3_M_eq_M(&U[d][idx], &CM[d]);
        } else {
            QMP_msghandle_t mh = QMP_declare_send_to(mm, s_node, 0);
            QMP_start(mh);
            QMP_wait(mh);
            QMP_free_msghandle(mh);
        }
    }
    QMP_free_msgmem(mm);
    fclose(f);

    if (cksum29 != hdrbuf.cksum[0] || cksum31 != hdrbuf.cksum[1]) {
        status = "checksum mismatch";
        goto clearerr_1;
    }

    stat_b      = broadcast_string(L, status);
    char *d = broadcast_string(L, datetmp);
    if (NULL != datetmp1) *datetmp1 = d;

    if (NULL != coord) qlua_free(L, coord);
    if (NULL != datetmp) qlua_free(L, datetmp);
    if (NULL != site_buf) qlua_free(L, site_buf);
    if (NULL != CM) qlua_free(L, CM);
    return 0;


clearerr_1:
clearerr_0:
    stat_b = broadcast_string(L, status);
    
    if (NULL != datetmp) qlua_free(L, datetmp);
    if (NULL != coord) qlua_free(L, coord);
    if (NULL != site_buf) qlua_free(L, site_buf);
    if (NULL != CM) qlua_free(L, CM);
    
    return luaL_error(L, status);
}

static int
read_milc_gauge_slave(lua_State *L,
                 mLattice *S,
                 char **datetmp1,
                 QLA_D3_ColorMatrix **U,
                 const char *name)
{
    char *status = NULL,
         *stat_b = NULL;

    long long volume;
    long long site;
    int i;
    QLA_D3_ColorMatrix *CM = qlua_malloc(L, S->rank * sizeof (QLA_D3_ColorMatrix));
    int *coord = qlua_malloc(L, S->rank * sizeof (int));
    QMP_msgmem_t mm = QMP_declare_msgmem(&CM[0], S->rank * sizeof (CM[0]));

    /* get gauge element for this node */
    for (volume = 1, i = 0; i < S->rank; i++)
        volume *= S->dim[i];

    /* get all data from node qlua_master_node */
    for (site = 0; site < volume; site++) {
        int s_node;
        site2coord(coord, site, S->rank, S->dim);
        s_node = QDP_node_number_L(S->lat, coord);
        if (s_node == QDP_this_node) {
            QMP_msghandle_t mh = QMP_declare_receive_from(mm, qlua_master_node, 0);
            QMP_start(mh);
            QMP_wait(mh);
            QMP_free_msghandle(mh);

            int idx = QDP_index_L(S->lat, coord);
            int d;
            
            for (d = 0; d < S->rank; d++)
                QLA_D3_M_eq_M(&U[d][idx], &CM[d]);
        }
    }
    if (NULL != mm) QMP_free_msgmem(mm);
    if (NULL != coord) qlua_free(L, coord);
    if (NULL != CM) qlua_free(L, CM);

    stat_b = broadcast_string(L, status);
    if (NULL != stat_b)
        return luaL_error(L, stat_b);

    char *d = broadcast_string(L, NULL);
    if (NULL != datetmp1)
        *datetmp1 = d;

    return 0;
}



static int
q_milc_gauge_read(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);
    const char *name = luaL_checkstring(L, 2);
    char *datetmp = NULL;
    int status;
    
    lua_createtable(L, S->rank, 0);
    QDP_D3_ColorMatrix **M = qlua_malloc(L, S->rank * sizeof (QDP_D3_ColorMatrix *));
    QLA_D3_ColorMatrix **U = qlua_malloc(L, S->rank * sizeof (QLA_D3_ColorMatrix *));
    if (NULL == M || NULL == U) {
        return luaL_error(L, "not enough memory");
    }
    CALL_QDP(L);
    for (int i = 0; i < S->rank; i++) {
        M[i] = qlua_newLatColMat3(L, 1, 3)->ptr;
        U[i] = QDP_D3_expose_M(M[i]);
        lua_rawseti(L, -2, i + 1);
    }
    
    if (QDP_this_node == qlua_master_node) {
        status = read_milc_gauge_master(L, S, &datetmp, U, name);
    } else {
        status = read_milc_gauge_slave(L, S, &datetmp, U, name);
    }
    CALL_QDP(L);
    for (int i = 0; i < S->rank; i++) {
        QDP_D3_reset_M(M[i]);
    }

    lua_newtable(L);
    lua_pushstring(L, datetmp);
    lua_setfield(L, -2, "datetmp");

    if (NULL != U) qlua_free(L, U);
    if (NULL != M) qlua_free(L, M);
    if (NULL != datetmp) qlua_free(L, datetmp);
    
    return 2;
}

static const struct luaL_Reg fMILC[] = {
    { "read_gauge",      q_milc_gauge_read },
    { NULL,              NULL         }
};

int
init_milc_io(lua_State *L)
{
    lua_getglobal(L, qcdlib);
    lua_newtable(L);
    luaL_register(L, NULL, fMILC);
    lua_setfield(L, -2, milc_io);
    lua_pop(L, 1);

    return 0;
}
#else /* USE_Nc3 */
int
init_milc_io(lua_State *L)
{
    return 0;
}
#endif /* USE_Nc3 == 0 */

void
fini_milc_io(void)
{
}
