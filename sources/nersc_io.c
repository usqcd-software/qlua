/* Only Nc=3 reader is supported */

#include "qlua.h"                                                    /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "latcolmat.h"                                               /* DEPS */
#include "qmp.h"
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

static const char nersc_io[] = "nersc";

enum {
    nerscERROR,  /* B: +size ++errstr[size] -- any error from the master */
    nerscOK,     /* B: +any -- file read successfully */
    nerscKEY,    /* B: +size ++(key[size]) -- header elem key */
    nerscVALUE,  /* B: +size ++(value[size]) -- header elem value */
    nerscEOH,    /* B: +any -- end of the header */
    nerscSTART,  /* B: +any -- start of data */
    nerscDATA    /* B: + 1 site gauge field -- gauge fields for one site */
};

typedef struct {
    int code;
    int arg;
} NERSC_Command;

static int
nersc_error(lua_State *L, const char *errmsg)
{
    return luaL_error(L, "L:read_NERSC_gauge() error: %s", errmsg);
}

static int
nersc_mpi_error(lua_State *L, QMP_status_t status)
{
    return nersc_error(L, QMP_error_string(status));
}

static int
nersc_master_error(lua_State *L, const char *msg)
{
    NERSC_Command cmd;
    QMP_status_t status;

    cmd.code = nerscERROR;
    cmd.arg = strlen(msg) + 1;
    status = QMP_broadcast(&cmd, sizeof (cmd));
    if (status != QMP_SUCCESS)
        return nersc_mpi_error(L, status);
    status = QMP_broadcast((char *)msg, cmd.arg);
    if (status != QMP_SUCCESS)
        return nersc_mpi_error(L, status);

    return nersc_error(L, msg);
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

static int
nersc_master_cmd(lua_State *L, int code, int arg)
{
    NERSC_Command cmd;
    QMP_status_t status;

    cmd.code = code;
    cmd.arg = arg;
    status = QMP_broadcast(&cmd, sizeof (cmd));
    if (status != QMP_SUCCESS)
        return nersc_mpi_error(L, status);
    return 0;
}

static int
nersc_master_bcast(lua_State *L, int code, int len, void *payload)
{
    QMP_status_t status;

    nersc_master_cmd(L, code, len);
    status = QMP_broadcast(payload, len);
    if (status != QMP_SUCCESS)
        return nersc_mpi_error(L, status);

    return 0;
}

typedef struct {
    const char *name;
    int value;
} NERSC_Value;

static int
decode_hdr(lua_State *L, const char *value, int old, int expected,
           const NERSC_Value *t, const char *msg)
{
    int i;

    if (old != expected)
        return nersc_master_error(L, msg);

    for (i = 0; t[i].name; i++) {
        if (strcmp(t[i].name, value) == 0)
            return t[i].value;
    }
    return nersc_master_error(L, msg);
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

    if (nd != S->rank)
        nersc_master_error(L, "internal error (read_3x3)");

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

    if (nd != S->rank)
        nersc_master_error(L, "internal error (read_3x2)");

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

static const char *ukey = "unitarity";

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
        /* nt4D_4x4 */
    };
    static const NERSC_Value nFMTs[] = {
        {"4D_SU3_GAUGE_3x3", nt4D_3x3 },
        {"4D_SU3_GAUGE",     nt4D_3x2 },
        /* {"4D_SU4_GAUGE",     nt4D_4x4 }, */
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
    int  coord[S->rank];
    long long volume;
    long long site;
    int s_node;
    char dim_ok[S->rank];
    int i;
    int f_format = ntNONE;
    int f_fp = 0;
    int f_cs_p = 0;
    uint32_t f_checksum = 0;
    uint32_t d_checksum = 0;
    GaugeReader read_matrix;
    RealReader read_real;
    int site_size;
    int big_endian;
    double uni_eps;
    QLA_D3_ColorMatrix mone;
    QLA_D_Complex cone;
    double max_eps = 0.0;

    QLA_c_eq_r_plus_ir(cone, 1.0, 0.0);
    QLA_D3_M_eq_c(&mone, &cone);

    if (f == 0)
        return nersc_master_error(L, "file open error");

    for (i = 0; i < S->rank; i++)
        dim_ok[i] = 0;

    /* parse the header and distribute it across the system */
    if ((nersc_gethdr(f, buffer, &key, &value) != 1) ||
        (strcmp(key, "BEGIN_HEADER") != 0))
        return nersc_master_error(L, "missing header");

    for (;;) {
        switch (nersc_gethdr(f, buffer, &key, &value)) {
        case 1:
            if (strcmp(key, "END_HEADER") != 0)
                return nersc_master_error(L, "missing end of header");
            nersc_master_cmd(L, nerscEOH, 0);
            goto eoh;
        case 2:
            lua_pushstring(L, value);
            lua_setfield(L, -2, key);

            nersc_master_bcast(L, nerscVALUE, strlen(value) + 1, value);
            nersc_master_bcast(L, nerscKEY, strlen(key) + 1, key);
            
            if (strcmp(key, "DATATYPE") == 0) {
                f_format = decode_hdr(L, value, f_format, ntNONE, nFMTs,
                                      "bad or conflicting DATATYPE");
            } else if (strcmp(key, "FLOATING_POINT") == 0) {
                f_fp = decode_hdr(L, value, f_fp, 0, nFPs,
                                  "bad or conflicting FLOATING_POINT");
            } else if (strcmp(key, "CHECKSUM") == 0) {
                if (f_cs_p)
                    return nersc_master_error(L, "multiple CHECKSUMs");
                if (sscanf(value, "%x", &f_checksum) != 1)
                    return nersc_master_error(L, "illformed CHECKSUM");
                f_cs_p = 1;
                lua_pushnumber(L, f_checksum);
                lua_setfield(L, -2, key);
            } else if (sscanf(key, "DIMENSION_%d", &i) == 1) {
                int di;
                if ((i < 1) || (i > S->rank))
                    return nersc_master_error(L, "DIMENSION out of range");
                if ((sscanf(value, "%d", &di) != 1) ||
                    (S->dim[i - 1] != di))
                    return nersc_master_error(L, "DIMENSION mismatch");
                dim_ok[i - 1] = 1;
                lua_pushnumber(L, di);
                lua_setfield(L, -2, key);
            } else if (sscanf(key, "BOUNDARY_%d", &i) == 1) {
                if ((i < 1) || (i > S->rank))
                    return nersc_master_error(L, "BOUNDARY out of range");
                if (strcmp(value, "PERIODIC") != 0)
                    return nersc_master_error(L, "bad BOUNDARY value");
            } else if ((strcmp(key, "PLAQUETTE") == 0) ||
                       (strcmp(key, "LINK_TRACE") == 0)) {
                double v;
                if (sscanf(value, "%lg", &v) != 1)
                    return nersc_master_error(L, "unexpected value");
                lua_pushnumber(L, v);
                lua_setfield(L, -2, key);
            }
            break;
        default:
            return nersc_master_error(L, "illformed header");
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
#if 0 /* XXX handling SU(4) */
    case nt4D_4x4:
        if (QDP_Nc != 4)
            return nersc_master_error(L, "Unsupported Nc");
        read_matrix = read_NxN;
        site_size = S->rank * 4 * 4 * 2;
        break;
#endif /* XXX handling SU(4) */
    default:
        return nersc_master_error(L, "unsupported data format");
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
        return nersc_master_error(L, "bad floating point size");
    }
    for (i = 0; i < S->rank; i++) {
        if (!dim_ok[i])
            return nersc_master_error(L, "missing DIMENSION spec");
    }
    if (f_cs_p == 0)
        return nersc_master_error(L, "missing CHECKSUM");

    /* Read the data and send it to the target host */
    for (volume = 1, i = 0; i < S->rank; i++)
        volume *= S->dim[i];
    char site_buf[site_size];
    QLA_D3_ColorMatrix CM[S->rank];

    /* find out our endianess */
    {
        union {
            uint64_t ll;
            unsigned char c[sizeof (uint64_t)];
        } b;
        uint64_t v;
        for (v = 1, i = 0; i < sizeof (uint64_t); i++)
            v = (v << CHAR_BIT) + i + 1;
        b.ll = v;
        if (b.c[0] == 1)
            big_endian = 1;
        else if (b.c[0] == i)
            big_endian = 0;
        else
            return nersc_master_error(L, "Unexpected host endianness");
    }

    /* start of data message */
    nersc_master_cmd(L, nerscSTART, 0);

    /* Go through all sites */
    for (site = 0; site < volume; site++) {
        if (fread(site_buf, site_size, 1, f) != 1)
            return nersc_master_error(L, "file read error");
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
                return nersc_master_error(L, "internal error");
            }
        }
        /* collect the checksum */
        for (i = 0; i < site_size; i += sizeof (uint32_t))
            d_checksum += *(uint32_t *)(site_buf + i);
        /* convert to the ColorMatrix */
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
                        if ((er > uni_eps) || (ei > uni_eps))
                            return nersc_master_error(L, "unitarty violation");
                        if (er > max_eps)
                            max_eps = er;
                        if (ei > max_eps)
                            max_eps = ei;
                    }
                }
            }
        }
        /* send CM to everyone
         *  -- this makes all slaves receive previous errors as well.
         */
        nersc_master_bcast(L, nerscDATA, sizeof (CM), CM);
        /* place ColorMatrix in U if it belongs to this site */
        site2coord(coord, site, S->rank, S->dim);
        s_node = QDP_node_number_L(S->lat, coord);
        if (s_node == QDP_this_node) {
            int idx = QDP_index_L(S->lat, coord);
            int d;
            
            for (d = 0; d < S->rank; d++)
                QLA_D3_M_eq_M(&U[d][idx], &CM[d]);
        }
    }

    fclose(f);

    /* check the checksum, broadcast nerscOK or nerscERROR */
    if (f_checksum != d_checksum)
        return nersc_master_error(L, "checksum mismatch");

    lua_pushnumber(L, max_eps);
    lua_setfield(L, -2, ukey);
    nersc_master_bcast(L, nerscVALUE, sizeof (double), &max_eps);
    nersc_master_cmd(L, nerscOK, 0);

    return 2;
}

static int
nersc_slave_error(lua_State *L, int len)
{
    char *msg = qlua_malloc(L, len);
    QMP_status_t status;

    status = QMP_broadcast(msg, len);
    if (status != QMP_SUCCESS)
        return nersc_mpi_error(L, status);

    return nersc_error(L, msg);
}

static int
nersc_slave_cmd(lua_State *L, int *code, int *arg)
{
    NERSC_Command cmd;
    QMP_status_t status = QMP_broadcast(&cmd, sizeof (cmd));

    if (status != QMP_SUCCESS)
        return nersc_mpi_error(L, status);

    if (code == nerscERROR)
        return nersc_slave_error(L, *arg);

    *code = cmd.code;
    *arg = cmd.arg;

    return cmd.code;
}

static int
nersc_slave_bcast(lua_State *L, int len, void *buffer)
{
    QMP_status_t status = QMP_broadcast(buffer, len);
    
    if (status != QMP_SUCCESS)
        return nersc_mpi_error(L, status);

    return 0;
}

static int
nersc_read_slave(lua_State *L,
                 mLattice *S,
                 QLA_D3_ColorMatrix **U,
                 const char *name)
{
    int code = 0;
    int arg = 0;
    char key[NERSC_BUFSIZE];
    char value[NERSC_BUFSIZE];
    long long volume;
    long long site;
    int coord[S->rank];
    int s_node;
    int i;
    QLA_D3_ColorMatrix CM[S->rank];
    double max_eps;

    /* get header elements */
    for (;;) {
        switch (nersc_slave_cmd(L, &code, &arg)) {
        case nerscEOH:
            goto eoh;
        case nerscVALUE:
            nersc_slave_bcast(L, arg, value);
            lua_pushstring(L, value);
            if (nersc_slave_cmd(L, &code, &arg) != nerscKEY)
                return luaL_error(L, "internal error (key %d)", code);
            nersc_slave_bcast(L, arg, key);
            lua_setfield(L, -2, key);
            if ((strcmp(key, "PLAQUETTE") == 0) ||
                (strcmp(key, "LINK_TRACE") == 0)) {
                double v = 0;
                sscanf(value, "%lg", &v); /* errors handled by the master */
                lua_pushnumber(L, v);
                lua_setfield(L, -2, key);
            } else if (strcmp(key, "CHECKSUM") == 0) {
                uint32_t s;
                sscanf(value, "%x", &s); /* errors handled by the master */
                lua_pushnumber(L, s);
                lua_setfield(L, -2, key);
            } else if (sscanf(key, "DIMENSION_%d", &i) == 1) {
                int d;
                sscanf(value, "%d", &d); /* errors handled by the master */
                lua_pushnumber(L, d);
                lua_setfield(L, -2, key);
            }
            break;
        default:
            return luaL_error(L, "internal error (header %d)", code);
        }
    }
eoh:
    switch (nersc_slave_cmd(L, &code, &arg)) {
    case nerscSTART:
        break;
    default:
        return luaL_error(L, "internal error (START %d)", code);
    }
    /* get gauge element for this node */
    for (volume = 1, i = 0; i < S->rank; i++)
        volume *= S->dim[i];

    for (site = 0; site < volume; site++) {
        switch (nersc_slave_cmd(L, &code, &arg)) {
        case nerscDATA:
            break;
        default:
            return luaL_error(L, "internal error (DATA %d)", code);
        }
        nersc_slave_bcast(L, sizeof (CM), CM);
        site2coord(coord, site, S->rank, S->dim);
        s_node = QDP_node_number_L(S->lat, coord);
        if (s_node == QDP_this_node) {
            int idx = QDP_index_L(S->lat, coord);
            int d;
            
            for (d = 0; d < S->rank; d++)
                QLA_D3_M_eq_M(&U[d][idx], &CM[d]);
        }
    }

    /* get the maximal unitary violation */
    switch (nersc_slave_cmd(L, &code, &arg)) {
    case nerscVALUE:
        if (arg != sizeof (double))
            luaL_error(L, "internal error (unitarity size)");
        break;
    default:
        return luaL_error(L, "internal error (UNITARITY %d)", code);
    }
    nersc_slave_bcast(L, sizeof (double), &max_eps);
    
    lua_pushnumber(L, max_eps);
    lua_setfield(L, -2, ukey);

    /* wait for the final nerscOK */
    switch (nersc_slave_cmd(L, &code, &arg)) {
    case nerscOK:
        break;
    default:
        return luaL_error(L, "internal error (END %d)", code);
    }

    return 2;
}

static int
q_nersc_read(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);
    const char *name = luaL_checkstring(L, 2);
    QDP_D3_ColorMatrix *M[S->rank];
    QLA_D3_ColorMatrix *U[S->rank];
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
    if (qlua_primary_node) {
        status = nersc_read_master(L, S, U, name);
    } else {
        status = nersc_read_slave(L, S, U, name);
    }
    CALL_QDP(L);
    for (i = 0; i < S->rank; i++) {
        QDP_D3_reset_M(M[i]);
    }
    
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

int
fini_nersc_io(lua_State *L)
{
    return 0;
}
