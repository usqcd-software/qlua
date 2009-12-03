#include <qlua.h>                                                    /* DEPS */
#include <qdpcc_io.h>                                                /* DEPS */
#include <qio_utils.h>                                               /* DEPS */
#include <lattice.h>                                                 /* DEPS */
#include <latdirferm.h>                                              /* DEPS */
#include <latdirprop.h>                                              /* DEPS */
#include <qio.h>

static const char qdpcc_io[] = "qdpcc";

#if QDP_Ns == 4
typedef QLA_F_Complex qdpcc_Ft[QDP_Ns][QDP_Ns][QDP_Nc][QDP_Nc];
typedef QLA_D_Complex qdpcc_Dt[QDP_Ns][QDP_Ns][QDP_Nc][QDP_Nc];
#else
/* define respective Chroma type here */
#endif

typedef struct {
    lua_State *L;
    QLA_DiracPropagator *P;
} QDPCCArgs;

static void
float_P_put(char *buf, size_t index, int count, void *arg0)
{
    qdpcc_Ft *src = (qdpcc_Ft *)buf;
    QDPCCArgs *arg = (QDPCCArgs *)arg0;
    int ic, jc, is, js;
    QLA_DiracPropagator *P = arg->P + index;

    if (count != 1)
        luaL_error(arg->L, "qcd.qdpcc.read_prop() count != 1");

    for (is = 0; is < QDP_Ns; is++) {
        for (js = 0; js < QDP_Ns; js++) {
            for (ic = 0; ic < QDP_Nc; ic++) {
                for (jc = 0; jc < QDP_Nc; jc++) {
                    QLA_F_Complex *z = &(*src)[is][js][ic][jc];
                    QLA_c_eq_r_plus_ir(QLA_elem_P(*P, ic, is, jc, js),
                                       QLA_real(*z), QLA_imag(*z));
                               
                }
            }
        }
    }
}

static void
double_P_put(char *buf, size_t index, int count, void *arg0)
{
    qdpcc_Dt *src = (qdpcc_Dt *)buf;
    QDPCCArgs *arg = (QDPCCArgs *)arg0;
    int ic, jc, is, js;
    QLA_DiracPropagator *P = arg->P + index;

    if (count != 1)
        luaL_error(arg->L, "qcd.qdpcc.read_prop() count != 1");

    for (is = 0; is < QDP_Ns; is++) {
        for (js = 0; js < QDP_Ns; js++) {
            for (ic = 0; ic < QDP_Nc; ic++) {
                for (jc = 0; jc < QDP_Nc; jc++) {
                    QLA_D_Complex *z = &(*src)[is][js][ic][jc];
                    QLA_c_eq_r_plus_ir(QLA_elem_P(*P, ic, is, jc, js),
                                       QLA_real(*z), QLA_imag(*z));
                               
                }
            }
        }
    }
}

static int
qdpcc_read_P(lua_State *L)
{
    typedef void Putter(char *buf, size_t index, int count, void *arg);
    static const struct {
        char    fmt;
        Putter *putter;
        int     word_size;
        size_t  datum_size;
    } read_P_fmt[] = {
        {'F', float_P_put,  sizeof (QLA_F_Real),sizeof (QLA_F_DiracPropagator)},
        {'D', double_P_put, sizeof (QLA_D_Real),sizeof (QLA_D_DiracPropagator)},
        {0, NULL, 0, 0}
    };
    const char *fname = luaL_checkstring(L, 1);
    mLatDirProp *P = qlua_newLatDirProp(L);
    QIO_Reader *reader = 0;
    QIO_String *file_xml = 0;
    QIO_String *record_xml = 0;
    QIO_RecordInfo *record_info = 0;
    Putter *putter = 0;
    size_t datum_size = 0;
    int word_size = 0;
    QDPCCArgs data;
    int i, fmt;

    CALL_QDP(L);
    file_xml = QIO_string_create();
    record_xml = QIO_string_create();
    reader = qlua_qio_std_reader(fname, file_xml);
    if (reader == 0)
        return luaL_error(L, "qcd.qdpcc.read_prop() open error");
    record_info = QIO_create_record_info(0, NULL, NULL, 0, "", "", 0, 0, 0, 0);
    if (QIO_read_record_info(reader, record_info, record_xml))
        return luaL_error(L, "qcd.qdpcc.read_prop() read info error");
    fmt = *QIO_get_precision(record_info);
    for (i = 0; read_P_fmt[i].putter; i++) {
        if (fmt == read_P_fmt[i].fmt) {
            putter = read_P_fmt[i].putter;
            datum_size = read_P_fmt[i].datum_size;
            word_size = read_P_fmt[i].word_size;
            break;
        }
    }
    if (putter == 0)
        luaL_error(L, "qcd.qdpcc.read_prop(): Unknown precision");
    data.P = QDP_expose_P(P->ptr);
    data.L = L;
    if (QIO_read(reader, record_info, record_xml,
                 putter, datum_size, word_size, &data))
        luaL_error(L, "qcd.qdpcc.read_prop() read data error");
    QDP_reset_P(P->ptr);
    QIO_destroy_record_info(record_info);
    QIO_close_read(reader);
    lua_pushstring(L, QIO_string_ptr(file_xml));
    lua_pushstring(L, QIO_string_ptr(record_xml));
    QIO_string_destroy(file_xml);
    QIO_string_destroy(record_xml);

    return 3;
}

static void
float_P_get(char *buf, size_t index, int count, void *arg0)
{
    QDPCCArgs *arg = (QDPCCArgs *)arg0;
    QLA_DiracPropagator *P = arg->P + index;
    qdpcc_Ft *dst = (qdpcc_Ft *)buf;
    int ic, jc, is, js;

    if (count != 1)
        luaL_error(arg->L, "qcd.qdpcc.write_prop() count != 1");

    for (is = 0; is < QDP_Ns; is++) {
        for (js = 0; js < QDP_Ns; js++) {
            for (ic = 0; ic < QDP_Nc; ic++) {
                for (jc = 0; jc < QDP_Nc; jc++) {
                    QLA_Complex *z = &QLA_elem_P(*P, ic, is, jc, is);
                    QLA_c_eq_r_plus_ir((*dst)[is][js][ic][jc],
                                       QLA_real(*z), QLA_imag(*z));
                }
            }
        }
    }
}

static void
double_P_get(char *buf, size_t index, int count, void *arg0)
{
    QDPCCArgs *arg = (QDPCCArgs *)arg0;
    QLA_DiracPropagator *P = arg->P + index;
    qdpcc_Dt *dst = (qdpcc_Dt *)buf;
    int ic, jc, is, js;

    if (count != 1)
        luaL_error(arg->L, "qcd.qdpcc.write_prop() count != 1");

    for (is = 0; is < QDP_Ns; is++) {
        for (js = 0; js < QDP_Ns; js++) {
            for (ic = 0; ic < QDP_Nc; ic++) {
                for (jc = 0; jc < QDP_Nc; jc++) {
                    QLA_Complex *z = &QLA_elem_P(*P, ic, is, jc, is);
                    QLA_c_eq_r_plus_ir((*dst)[is][js][ic][jc],
                                       QLA_real(*z), QLA_imag(*z));
                }
            }
        }
    }
}

static int
qdpcc_write_P(lua_State *L)
{
    typedef void Getter(char *buf, size_t index, int count, void *arg);
#if QDP_Nc == 3
    static const struct {
        char *fmt;
        Getter *getter;
        char *datatype;
        int word_size;
        int datum_size;
    } write_P_fmt[] = {
        {"F", float_P_get, "QDP_F3_DiracPropagator",
         sizeof (QLA_F_Real), sizeof (QLA_F_DiracPropagator)},
        {"D", double_P_get, "QDP_D3_DiracPropagator",
         sizeof (QLA_D_Real), sizeof (QLA_D_DiracPropagator)},
        { 0, NULL, NULL,
          0, 0}
    };
#else
    /* define write_P_fmt here */
#endif
    const char *fmt = luaL_checkstring(L, 1);
    const char *fname = luaL_checkstring(L, 2);
    const char *file_xml = luaL_checkstring(L, 3);
    mLatDirProp *P = qlua_checkLatDirProp(L, 4);
    const char *record_xml = luaL_checkstring(L, 5);
    QIO_Writer *writer = NULL;
    Getter *getter = NULL;
    char *prec = NULL;
    int typesize = 0;
    int datacount = 0;
    int datum_size = 0;
    int word_size = 0;
    char *datatype = NULL;
    QIO_String *f_xml = NULL;
    QIO_String *r_xml = NULL;
    QIO_RecordInfo *record_info = NULL;
    QDPCCArgs data;
    int i;

    CALL_QDP(L);
    for (i = 0; write_P_fmt[i].getter; i++) {
        if (fmt[0] == write_P_fmt[i].fmt[0]) {
            prec = write_P_fmt[i].fmt;
            getter = write_P_fmt[i].getter;
            datatype = write_P_fmt[i].datatype;
            datum_size = write_P_fmt[i].datum_size;
            word_size = write_P_fmt[i].word_size;
            datacount = 1;
            typesize = datum_size;
            break;
        }
    }
    if (getter == 0)
        return luaL_error(L, "qcd.qdpcc.write_prop(): Unknown precission");
    f_xml = QIO_string_create();
    QIO_string_set(f_xml, file_xml);
    r_xml = QIO_string_create();
    QIO_string_set(r_xml, record_xml);
    writer = qlua_qio_std_writer(fname, f_xml);
    if (writer == 0)
        return luaL_error(L, "qcd.qdpcc.write_prop(): open error");
    record_info = QIO_create_record_info(QIO_FIELD, NULL, NULL, qRank,
                                         datatype, prec, QDP_Nc, QDP_Ns,
                                         typesize, datacount);
    data.P = QDP_expose_P(P->ptr);
    data.L = L;
    if (QIO_write(writer, record_info, r_xml, getter,
                  datum_size, word_size, &data))
        return luaL_error(L, "qcd.qdpcc.write_prop(): record writing error");
        
    if (QIO_close_write(writer))
        return luaL_error(L, "qcd.qdpcc.write_prop(): closing error");

    QDP_reset_P(P->ptr);
    QIO_destroy_record_info(record_info);
    QIO_string_destroy(f_xml);
    QIO_string_destroy(r_xml);

    return 0;
}

static struct luaL_Reg fQDPCCio[] = {
    { "read_prop",     qdpcc_read_P },
    { "write_prop",    qdpcc_write_P },
    { NULL,       NULL }
};

int
init_qdpcc_io(lua_State *L)
{
    int i;

    lua_getglobal(L, qcdlib);
    lua_newtable(L);
    for (i = 0; fQDPCCio[i].name; i++) {
        lua_pushcfunction(L, fQDPCCio[i].func);
        lua_setfield(L, -2, fQDPCCio[i].name);
    }
    lua_setfield(L, -2, qdpcc_io);
    lua_pop(L, 1);
    return 0;
}

int
fini_qdpcc_io(lua_State *L)
{
    return 0;
}
