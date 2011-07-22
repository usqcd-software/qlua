/* Only Nc=3 reader is supported */
#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "latcolmat.h"                                               /* DEPS */
#include "qmp.h"
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#if USE_Nc3
static const char nersc_io[] = "nersc";

static int
nersc_error(lua_State *L, const char *errmsg)
{
    return luaL_error(L, "L:read_NERSC_gauge() error: %s", errmsg);
}

static char *
nersc_hdrnorm(char *buffer)
{
    char *head, *tail;

    head = buffer;
    tail = buffer + strlen(buffer) - 1;
    
    while (*head == ' ')
        ++head;
    while (*tail == ' ')
        *tail-- = 0;

    return head;
}

#define NERSC_BUFSIZE 1024

static int
nersc_gethdr(FILE *f, char *buffer, char **key, char **value)
{
    char *eq;

    if (fgets(buffer, NERSC_BUFSIZE, f) == 0)
        return 0;
    buffer[NERSC_BUFSIZE - 1] = 0;
    eq = buffer + strlen(buffer) - 1;
    if (*eq != '\n')
        return 0;
    *eq = 0;

    eq = strchr(buffer, '=');
    if (eq == NULL) {
        *key = nersc_hdrnorm(buffer);
        *value = NULL;
        return 1;
    } else {
        *eq++ = 0;
        *key = nersc_hdrnorm(buffer);
        *value = nersc_hdrnorm(eq);
        return 2;
    }
}

typedef struct {
    const char *name;
    int value;
} NERSC_Value;

static int
decode_hdr(lua_State *L, const char *value, int old, int expected,
           const NERSC_Value *t, char *msg, char **status)
{
    int i;

    if (*status != NULL)
        return old;

    if (old != expected) {
        *status = msg;
        return old;
    }

    for (i = 0; t[i].name; i++) {
        if (strcmp(t[i].name, value) == 0)
            return t[i].value;
    }
    *status = msg;
    return old;
}

typedef double (*RealReader)(char *data, int idx);
static double
read_float(char *data, int idx)
{
    return *(float *)(data + idx * sizeof (float));
}

static double
read_double(char *data, int idx)
{
    return *(double *)(data + idx * sizeof (double));
}

typedef void (*GaugeReader)(lua_State *L,
                            mLattice *S,
                            QLA_D3_ColorMatrix *U, int nd,
                            char *buf, int buf_size,
                            RealReader read_real);

static void
read_3x3(lua_State *L,
         mLattice *S,
         QLA_D3_ColorMatrix *U, int nd,
         char *buf, int buf_size,
         RealReader read_real)
{
    int i, d, a, b;

    for (i = 0, d = 0; d < nd; d++) {
        for (a = 0; a < 3; a++) {
            for (b = 0; b < 3; b++, i += 2) {
                QLA_c_eq_r_plus_ir(QLA_D3_elem_M(U[d], a, b),
                                   read_real(buf, i),
                                   read_real(buf, i + 1));
            }
        }
    }
}

static void
read_3x2(lua_State *L,
         mLattice *S,
         QLA_D3_ColorMatrix *U, int nd,
         char *buf, int buf_size,
         RealReader read_real)
{
    int i, d, a, b;

    for (i = 0, d = 0; d < nd; d++, U++) {
        for (a = 0; a < 3 - 1; a++) {
            for (b = 0; b < 3; b++, i += 2) {
                QLA_c_eq_r_plus_ir(QLA_D3_elem_M(*U, a, b),
                                   read_real(buf, i),
                                   read_real(buf, i + 1));
            }
        }
        QLA_c_eq_ca_times_ca( QLA_D3_elem_M(*U,2,2), QLA_D3_elem_M(*U,0,0),
                              QLA_D3_elem_M(*U,1,1));
        QLA_c_meq_ca_times_ca(QLA_D3_elem_M(*U,2,2), QLA_D3_elem_M(*U,0,1),
                              QLA_D3_elem_M(*U,1,0));

        QLA_c_eq_ca_times_ca( QLA_D3_elem_M(*U,2,1), QLA_D3_elem_M(*U,0,2),
                              QLA_D3_elem_M(*U,1,0));
        QLA_c_meq_ca_times_ca(QLA_D3_elem_M(*U,2,1), QLA_D3_elem_M(*U,0,0),
                              QLA_D3_elem_M(*U,1,2));

        QLA_c_eq_ca_times_ca( QLA_D3_elem_M(*U,2,0), QLA_D3_elem_M(*U,0,1),
                              QLA_D3_elem_M(*U,1,2));
        QLA_c_meq_ca_times_ca(QLA_D3_elem_M(*U,2,0), QLA_D3_elem_M(*U,0,2),
                              QLA_D3_elem_M(*U,1,1));
    }
}

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
send_string(lua_State *L, int idx)
{
    const char *str = lua_tostring(L, idx);
    int len = strlen(str) + 1;

    QMP_broadcast(&len, sizeof (len));
    QMP_broadcast((void *)str, len);
}

static int
receive_string(lua_State *L, char **str)
{
    int len;
    *str = 0;
    QMP_broadcast(&len, sizeof (len));
    if (len == 0)
        return 0;
    *str = qlua_malloc(L, len);
    QMP_broadcast(*str, len);
    return 1;
}

static void
normalize_int(lua_State *L, int idx, const char *key, char *fmt)
{
    lua_getfield(L, idx, key);
    const char *value = lua_tostring(L, -1);
    int v;
    if (value == NULL)
        return;

    if (sscanf(value, fmt, &v) == 1) {
        lua_pop(L, 1);
        lua_pushnumber(L, v);
        lua_setfield(L, idx - 1, key);
        return;
    }
    lua_pop(L, 1);
}

static void
normalize_float(lua_State *L, int idx, const char *key)
{
    lua_getfield(L, idx, key);
    const char *value = lua_tostring(L, -1);
    double v;
    if (value == NULL)
        return;
    if (sscanf(value, "%lg", &v) == 1) {
        lua_pop(L, 1);
        lua_pushnumber(L, v);
        lua_setfield(L, idx - 1, key);
        return;
    }
    lua_pop(L, 1);
}

static const char *ukey = "unitarity";

static void
normalize_kv(lua_State *L, mLattice *S, int idx)
{
    int i;
    for (i = 1; i <= S->rank; i++) {
        char buf[128]; /* large enough for DIMENSION_%d */
        snprintf(buf, sizeof (buf) - 1, "DIMENSION_%d", i);
        normalize_int(L, idx, buf, "%d");
    }
    normalize_int(L, idx, "CHECKSUM", "%x");
    normalize_float(L, idx, "PLAQUETTE");
    normalize_float(L, idx, "LINK_TRACE");
    normalize_float(L, idx, ukey);
}

static int
nersc_read_master(lua_State *L,
                  mLattice *S,
                  QLA_D3_ColorMatrix **U,
                  const char *name)
{
    enum {
        ntNONE,
        nt4D_3x3,
        nt4D_3x2
    };
    static const NERSC_Value nFMTs[] = {
        {"4D_SU3_GAUGE_3x3", nt4D_3x3 },
        {"4D_SU3_GAUGE",     nt4D_3x2 },
        {NULL, ntNONE}
    };
    static const NERSC_Value nFPs[] = {
        {"IEEE32",    4},
        {"IEEE32BIG", 4},
        {"IEEE64BIG", 8},
        {NULL,        0}
    };

    FILE *f = fopen(name, "rb");
    char buffer[NERSC_BUFSIZE];
    char *key, *value;
    int  *coord = qlua_malloc(L, S->rank * sizeof (int));
    char *dim_ok = qlua_malloc(L, S->rank * sizeof (int));
    long long volume;
    long long site;
    int s_node;
    int i;
    int f_format = ntNONE;
    int f_fp = 0;
    int f_cs_p = 0;
    uint32_t f_checksum = 0;
    uint32_t d_checksum = 0;
    GaugeReader read_matrix = NULL;
    RealReader read_real = NULL;
    int site_size = 0;
    int big_endian;
    double uni_eps = 0;
    QLA_D3_ColorMatrix mone;
    QLA_D_Complex cone;
    double max_eps = 0.0;
    char *status = NULL;

    QLA_c_eq_r_plus_ir(cone, 1.0, 0.0);
    QLA_D3_M_eq_c(&mone, &cone);

    if (f == 0)
        status = "file open error";

    for (i = 0; i < S->rank; i++)
        dim_ok[i] = 0;

    /* parse the header */
    if ((status == NULL) && 
        ((nersc_gethdr(f, buffer, &key, &value) != 1) ||
         (strcmp(key, "BEGIN_HEADER") != 0)))
        status = "missing header";

    while (status == NULL) {
        switch (nersc_gethdr(f, buffer, &key, &value)) {
        case 1:
            if (strcmp(key, "END_HEADER") != 0)
                status = "missing end of header";
            goto eoh;
        case 2:
            lua_pushstring(L, value);
            lua_setfield(L, -2, key);
            
            if (strcmp(key, "DATATYPE") == 0) {
                f_format = decode_hdr(L, value, f_format, ntNONE, nFMTs,
                                      "bad or conflicting DATATYPE", &status);
            } else if (strcmp(key, "FLOATING_POINT") == 0) {
                f_fp = decode_hdr(L, value, f_fp, 0, nFPs,
                                  "bad or conflicting FLOATING_POINT", &status);
            } else if (strcmp(key, "CHECKSUM") == 0) {
                if (f_cs_p)
                    status = "multiple CHECKSUMs";
                if ((status == NULL) && 
                    (sscanf(value, "%x", &f_checksum) != 1))
                    status = "illformed CHECKSUM";
                f_cs_p = 1;
            } else if (sscanf(key, "DIMENSION_%d", &i) == 1) {
                int di;
                if ((i < 1) || (i > S->rank))
                    status = "DIMENSION out of range";
                if ((status == NULL) &&
                    ((sscanf(value, "%d", &di) != 1) ||
                     (S->dim[i - 1] != di)))
                    status = "DIMENSION mismatch";
                dim_ok[i - 1] = 1;
            } else if (sscanf(key, "BOUNDARY_%d", &i) == 1) {
                if ((i < 1) || (i > S->rank))
                    status = "BOUNDARY out of range";
                if ((status == NULL) &&
                    (strcmp(value, "PERIODIC") != 0))
                    status = "bad BOUNDARY value";
            } else if ((strcmp(key, "PLAQUETTE") == 0) ||
                       (strcmp(key, "LINK_TRACE") == 0)) {
                double v;
                if (sscanf(value, "%lg", &v) != 1)
                    status = "unexpected value";
            }
            break;
        default:
            status = "illformed header";
        }
    }
eoh:
    switch (f_format) {
    case nt4D_3x3:
        read_matrix = read_3x3;
        site_size = S->rank * 3 * 3 * 2;
        break;
    case nt4D_3x2:
        read_matrix = read_3x2;
        site_size = S->rank * 3 * (3 - 1) * 2;
        break;
    default:
        if (status == NULL)
            status = "unsupported data format";
    }
    switch (f_fp) {
    case 4:
        read_real = read_float;
        site_size *= 4;
        uni_eps = 1e-6;
        break;
    case 8:
        read_real = read_double;
        site_size *= 8;
        uni_eps = 1e-12;
        break;
    default:
        if (status == NULL)
            status = "bad floating point size";
    }
    for (i = 0; i < S->rank; i++) {
        if ((!dim_ok[i]) && (status == NULL)) {
            status = "missing DIMENSION spec";
        }
    }
    if ((f_cs_p == 0) && (status == NULL))
        status = "missing CHECKSUM";

    /* Read the data and send it to the target host */
    for (volume = 1, i = 0; i < S->rank; i++)
        volume *= S->dim[i];
    /* find out our endianess */
    {
        union {
            uint64_t ll;
            unsigned char c[sizeof (uint64_t)];
        } b;
        uint64_t v;
        big_endian = 1;
        for (v = 1, i = 0; i < sizeof (uint64_t); i++)
            v = (v << CHAR_BIT) + i + 1;
        b.ll = v;
        if (b.c[0] == 1)
            big_endian = 1;
        else if (b.c[0] == i)
            big_endian = 0;
        else if (status != NULL)
            status = "Unexpected host endianness";
    }
    /* read every site in order on the master
     * compute the checksum and send it to the target node
     */
    char *site_buf = qlua_malloc(L, site_size);
    QLA_D3_ColorMatrix *CM = qlua_malloc(L, S->rank * sizeof (QLA_D3_ColorMatrix));
    QMP_msgmem_t mm = QMP_declare_msgmem(&CM[0], S->rank * sizeof (CM[0]));

    /* Go through all sites */
    for (site = 0; site < volume; site++) {
        if ((status == NULL) && (fread(site_buf, site_size, 1, f) != 1))
            status = "file read error";

        /* swap bytes if necessary */
        if (big_endian == 0) {
            char *p;
            switch (f_fp) {
            case 4:
                for (i = 0, p = site_buf; i < site_size; i += 4, p += 4) {
                    char t;
                    t = p[0]; p[0] = p[3]; p[3] = t;
                    t = p[1]; p[1] = p[2]; p[2] = t;
                }
                break;
            case 8:
                for (i = 0, p = site_buf; i < site_size; i += 8, p += 8) {
                    char t;
                    t = p[0]; p[0] = p[7]; p[7] = t;
                    t = p[1]; p[1] = p[6]; p[6] = t;
                    t = p[2]; p[2] = p[5]; p[5] = t;
                    t = p[3]; p[3] = p[4]; p[4] = t;
                }
                break;
            default:
                if (status == NULL) 
                    status = "internal error: unsupported f_fp in endiannes conversion";
            }
        }
        /* collect the checksum */
        for (i = 0; i < site_size; i += sizeof (uint32_t))
            d_checksum += *(uint32_t *)(site_buf + i);
        /* convert to the ColorMatrix */
        if (read_matrix != NULL && read_real != NULL) 
            read_matrix(L, S, CM, S->rank, site_buf, site_size, read_real); 
        /* check unitarity */
        {
            int d, a, b;
            QLA_D3_ColorMatrix UxU;

            for (d = 0; d < S->rank; d++) {
                double er, ei;

                /* multiplcation order helps to detect reconstruction bugs */
                QLA_D3_M_eq_M_times_Ma(&UxU, &CM[d], &CM[d]);
                QLA_D3_M_meq_M(&UxU, &mone);
                for (a = 0; a < 3; a++) {
                    for (b = 0; b < 3; b++) {
                        QLA_Complex z;

                        QLA_D3_C_eq_elem_M(&z, &UxU, a, b);
                        er = fabs(QLA_real(z));
                        ei = fabs(QLA_imag(z));
                        if (((er > uni_eps) || (ei > uni_eps)) && (status == NULL))
                            status = "unitarty violation";
                        if (er > max_eps)
                            max_eps = er;
                        if (ei > max_eps)
                            max_eps = ei;
                    }
                }
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
    qlua_free(L, coord);
    qlua_free(L, dim_ok);
    qlua_free(L, site_buf);
    qlua_free(L, CM);
        
    fclose(f);

    /* check the checksum */
    if ((status == NULL) && (f_checksum != d_checksum))
        status = "checksum mismatch";

    /* broadcast the size of status to everyone */
    int status_len = (status != NULL)? (strlen(status) + 1): 0;
    QMP_sum_int(&status_len);
    /* if status is not NULL, broadcast to to everyone and call lua error */
    if (status_len != 0) {
        QMP_broadcast(status, status_len);
        return nersc_error(L, status);
    }

    /* convert max_eps to a string and store it (L, -2, ukey) */
    snprintf(buffer, sizeof (buffer) - 1, "%.10e", max_eps);
    lua_pushstring(L, buffer);
    lua_setfield(L, -2, ukey);
    /* ditribute the table across the machine */
    {
        lua_pushnil(L);
        while (lua_next(L, -2) != 0) {
            send_string(L, -2); /* key */
            send_string(L, -1); /* value */
            lua_pop(L, 1);
        }
        int zero = 0;
        QMP_broadcast(&zero, sizeof (zero));
    }
    /* normalize key/value pairs */
    normalize_kv(L, S, -1);
    return 2;
}

static int
nersc_read_slave(lua_State *L,
                 mLattice *S,
                 QLA_D3_ColorMatrix **U,
                 const char *name)
{
    long long volume;
    long long site;
    int i;
    QLA_D3_ColorMatrix *CM = qlua_malloc(L, S->rank * sizeof (QLA_D3_ColorMatrix));
    int *coord = qlua_malloc(L, S->rank * sizeof (int));
    QMP_msgmem_t mm = QMP_declare_msgmem(&CM[0], S->rank * sizeof (CM[0]));
    QMP_msghandle_t mh = QMP_declare_receive_from(mm, qlua_master_node, 0);

    /* get gauge element for this node */
    for (volume = 1, i = 0; i < S->rank; i++)
        volume *= S->dim[i];

    /* get all data from node qlua_master_node */
    for (site = 0; site < volume; site++) {
        int s_node;
        site2coord(coord, site, S->rank, S->dim);
        s_node = QDP_node_number_L(S->lat, coord);
        if (s_node == QDP_this_node) {
            QMP_start(mh);
            QMP_wait(mh);

            int idx = QDP_index_L(S->lat, coord);
            int d;
            
            for (d = 0; d < S->rank; d++)
                QLA_D3_M_eq_M(&U[d][idx], &CM[d]);
        }
    }
    QMP_free_msghandle(mh);
    QMP_free_msgmem(mm);
    qlua_free(L, coord);
    qlua_free(L, CM);

    /* get status size (broadcast) */
    int status_len = 0;
    QMP_sum_int(&status_len);
    /* if not zero, get the message, call lua_error */
    if (status_len != 0) {
        char *msg = qlua_malloc(L, status_len);
        QMP_broadcast(msg, status_len);
        return nersc_error(L, msg); /* leaks msg, but it's in the abort path only */
    }
        
    /* get keys and values */
    /* NB: This may cause a slave node to run out of memory out of sync with the master */
    for (;;) {
        char *key, *value;
        if (receive_string(L, &key) == 0)
            break;
        receive_string(L, &value);
        lua_pushstring(L, value);
        lua_setfield(L, -2, key);
        qlua_free(L, key);
        qlua_free(L, value);
    }
        
    /* normalize key/value table */
    normalize_kv(L, S, -1);
    return 2;
}


static int
q_nersc_read(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);
    const char *name = luaL_checkstring(L, 2);
    QDP_D3_ColorMatrix **M = qlua_malloc(L, S->rank * sizeof (QDP_D3_ColorMatrix *));
    QLA_D3_ColorMatrix **U = qlua_malloc(L, S->rank * sizeof (QLA_D3_ColorMatrix *));
    int status;
    int i;
    
    lua_createtable(L, S->rank, 0);
    CALL_QDP(L);
    for (i = 0; i < S->rank; i++) {
        M[i] = qlua_newLatColMat3(L, 1, 3)->ptr;
        U[i] = QDP_D3_expose_M(M[i]);
        lua_rawseti(L, -2, i + 1);
    }
    lua_newtable(L);
    if (QDP_this_node == qlua_master_node) {
        status = nersc_read_master(L, S, U, name);
    } else {
        status = nersc_read_slave(L, S, U, name);
    }
    CALL_QDP(L);
    for (i = 0; i < S->rank; i++) {
        QDP_D3_reset_M(M[i]);
    }
    qlua_free(L, U);
    qlua_free(L, M);
    
    return status;
}

static const struct luaL_Reg fNERSC[] = {
    { "read_gauge",      q_nersc_read },
    { NULL,              NULL         }
};

int
init_nersc_io(lua_State *L)
{
    lua_getglobal(L, qcdlib);
    lua_newtable(L);
    luaL_register(L, NULL, fNERSC);
    lua_setfield(L, -2, nersc_io);
    lua_pop(L, 1);

    return 0;
}
#else /* USE_Nc3 */
int
init_nersc_io(lua_State *L)
{
    return 0;
}
#endif /* USE_Nc3 == 0 */

int
fini_nersc_io(lua_State *L)
{
    return 0;
}
