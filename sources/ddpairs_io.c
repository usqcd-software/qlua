#include <qlua.h>                                                    /* DEPS */
#include <ddpairs_io.h>                                              /* DEPS */
#include <qio_utils.h>                                               /* DEPS */
#include <lattice.h>                                                 /* DEPS */
#include <latdirprop.h>                                              /* DEPS */
#include <qio.h>

static const char ddpairs_io[] = "ddpairs";

typedef QLA_F_Complex usqcd_Ft[QDP_Ns][QDP_Nc];
typedef QLA_D_Complex usqcd_Dt[QDP_Ns][QDP_Nc];

typedef struct {
    lua_State             *L;
    QLA_DiracPropagator   *P;
    int                    js;
    int                    jc;
} USQCDArgs;

/* reader */
static void
float_DD_put(char *buf, size_t index, int count, void *arg0)
{
    USQCDArgs *arg = (USQCDArgs *)arg0;
    usqcd_Ft *src = (usqcd_Ft *)buf;
    QLA_DiracPropagator *dst = &arg->P[index];
    int js = arg->js;
    int jc = arg->jc;
    int is, ic;
    
    if (count != 1)
        luaL_error(arg->L, "qcd.ddpairs.read(): count != 1");
    
    for (is = 0; is < QDP_Ns; is++) {
        for (ic = 0; ic < QDP_Nc; ic++) {
            QLA_F_Complex *z = &(*src)[is][ic];
            QLA_c_eq_r_plus_ir(QLA_elem_P(*dst, ic, is, jc, js),
                               QLA_real(*z), QLA_imag(*z));
        }
    }
}

static void
double_DD_put(char *buf, size_t index, int count, void *arg0)
{
    USQCDArgs *arg = (USQCDArgs *)arg0;
    usqcd_Dt *src = (usqcd_Dt *)buf;
    QLA_DiracPropagator *dst = &arg->P[index];
    int js = arg->js;
    int jc = arg->jc;
    int is, ic;
    
    if (count != 1)
        luaL_error(arg->L, "qcd.ddpairs.read(): count != 1");
    
    for (is = 0; is < QDP_Ns; is++) {
        for (ic = 0; ic < QDP_Nc; ic++) {
            QLA_D_Complex *z = &(*src)[is][ic];
            QLA_c_eq_r_plus_ir(QLA_elem_P(*dst, ic, is, jc, js),
                               QLA_real(*z), QLA_imag(*z));
        }
    }
}

static int
dd_read_comp(lua_State *L,
             int js, int jc,
             QLA_DiracPropagator *T,
             QIO_Reader *reader,
             QIO_RecordInfo *record_info,
             QIO_String *record_xml)
{
    typedef void Putter(char *buf, size_t index, int count, void *arg);
    static const struct {
        char    fmt;
        Putter *putter;
        int     word_size;
        size_t  datum_size;
    } read_DD_fmt[] = {
        {'F', float_DD_put,  sizeof (QLA_F_Real),sizeof (QLA_F_DiracFermion)},
        {'D', double_DD_put, sizeof (QLA_D_Real),sizeof (QLA_D_DiracFermion)},
        {0, NULL, 0, 0}
    };
    int word_size = 0;
    int datum_size = 0;
    Putter *putter = NULL;
    int i, prec;
    USQCDArgs data;

    if (QIO_read_record_info(reader, record_info, record_xml))
        return luaL_error(L, "qcd.ddpairs.read() read info error");
    prec = *QIO_get_precision(record_info);
    for (i = 0; read_DD_fmt[i].putter; i++) {
        if (prec == read_DD_fmt[i].fmt) {
            putter = read_DD_fmt[i].putter;
            word_size = read_DD_fmt[i].word_size;
            datum_size = read_DD_fmt[i].datum_size;
            break;
        }
    }
    if (putter == 0)
        return luaL_error(L, "qcd.ddpairs.read() unknown precision");
    data.P = T;
    data.js = js;
    data.jc = jc;
    data.L = L;
    if (QIO_read(reader, record_info, record_xml, putter,
                 datum_size, word_size, &data))
        return luaL_error(L, "qcd.ddpairs.read() read error");

    return 0;
}


static int
ddpairs_read(lua_State *L)
{
    const char *fname = luaL_checkstring(L, 1);
    mLatDirProp *S = qlua_newLatDirProp(L);
    mLatDirProp *P = qlua_newLatDirProp(L);
    QIO_String *file_xml = NULL;
    QIO_String *record_xml = NULL;
    QIO_RecordInfo *record_info = NULL;
    QIO_Reader *reader = NULL;
    int js, jc;
    QLA_DiracPropagator *xP;
    QLA_DiracPropagator *xS;
    
    lua_gc(L, LUA_GCCOLLECT, 0);
    file_xml = QIO_string_create();
    record_xml = QIO_string_create();
    record_info = QIO_create_record_info(0, NULL, NULL, 0, "", "", 0, 0, 0, 0);
    reader = qlua_qio_std_reader(fname, file_xml);
    if (reader == 0)
        return luaL_error(L, "qcd.ddpairs.read() open error");
    lua_pushstring(L, QIO_string_ptr(file_xml));
    QIO_string_destroy(file_xml);
    xP = QDP_expose_P(P->ptr);
    xS = QDP_expose_P(S->ptr);
    for (js = 0; js < QDP_Ns; js++) {
        for (jc = 0; jc < QDP_Nc; jc++) {
            dd_read_comp(L, js, jc, xS, reader, record_info, record_xml);
            dd_read_comp(L, js, jc, xP, reader, record_info, record_xml);
        }
    }
    QDP_reset_P(P->ptr);
    QDP_reset_P(S->ptr);
    QIO_destroy_record_info(record_info);
    QIO_string_destroy(record_xml);
    return 3;
}

/* writer */
static void
float_DD_get(char *buf, size_t index, int count, void *arg0)
{
    USQCDArgs *arg = (USQCDArgs *)arg0;
    QLA_DiracPropagator *P = arg->P + index;
    usqcd_Ft *dst = (usqcd_Ft *)buf;
    int jc = arg->jc;
    int js = arg->js;
    int ic, is;

    for (is = 0; is < QDP_Ns; is++) {
        for (ic = 0; ic < QDP_Nc; ic++) {
            QLA_Complex *z = &QLA_elem_P(*P, ic, is, jc, js);
            QLA_c_eq_r_plus_ir((*dst)[is][ic], QLA_real(*z), QLA_imag(*z));
        }
    }
}

static void
double_DD_get(char *buf, size_t index, int count, void *arg0)
{
    USQCDArgs *arg = (USQCDArgs *)arg0;
    QLA_DiracPropagator *P = arg->P + index;
    usqcd_Dt *dst = (usqcd_Dt *)buf;
    int jc = arg->jc;
    int js = arg->js;
    int ic, is;

    for (is = 0; is < QDP_Ns; is++) {
        for (ic = 0; ic < QDP_Nc; ic++) {
            QLA_Complex *z = &QLA_elem_P(*P, ic, is, jc, js);
            QLA_c_eq_r_plus_ir((*dst)[is][ic], QLA_real(*z), QLA_imag(*z));
        }
    }
}

static int
ddpairs_write(lua_State *L)
{
    typedef void Getter(char *buf, size_t index, int count, void *arg);
#if QDP_Nc == 3
    static const struct {
        char *fmt;
        Getter *getter;
        char *datatype;
        int word_size;
        int datum_size;
    } write_DD_fmt[] = {
        {"F", float_DD_get, "QDP_F3_DiracFermion",
         sizeof (QLA_F_Real), sizeof (QLA_F_DiracFermion)},
        {"D", double_DD_get, "QDP_D3_DiracFermion",
         sizeof (QLA_D_Real), sizeof (QLA_D_DiracFermion)},
        { 0, NULL, NULL,
          0, 0}
    };
#else
    /* define write_P_fmt here */
#endif
    const char *fmt = luaL_checkstring(L, 1);
    const char *fname = luaL_checkstring(L, 2);
    const char *file_xml = luaL_checkstring(L, 3);
    mLatDirProp *S = qlua_checkLatDirProp(L, 4);
    const char *source_xml = luaL_checkstring(L, 5);
    int time_slice = luaL_checkint(L, 6);
    mLatDirProp *P = qlua_checkLatDirProp(L, 7);
    const char *prop_xml = luaL_checkstring(L, 8);
    QIO_Writer *writer = NULL;
    Getter *getter = NULL;
    char *prec = NULL;
    int typesize = 0;
    int datacount = 0;
    int datum_size = 0;
    int word_size = 0;
    char *datatype = NULL;
    QIO_String *f_xml = NULL;
    QIO_String *s_xml = NULL;
    QIO_String *p_xml = NULL;
    USQCDArgs dataP;
    USQCDArgs dataS;
    int i, jc, js;
    int *lo = qlua_malloc(L, qRank * sizeof (int));
    int *hi = qlua_malloc(L, qRank * sizeof (int));

    lua_gc(L, LUA_GCCOLLECT, 0);
    for (i = 0; i < qRank; i++) {
        lo[i] = 0;
        hi[i] = qDim[i] - 1;
    }
    if ((time_slice >= 0) && (time_slice < qDim[qRank - 1]))
        lo[qRank - 1] = hi[qRank - 1] = time_slice;
    for (i = 0; write_DD_fmt[i].getter; i++) {
        if (fmt[0] == write_DD_fmt[i].fmt[0]) {
            prec = write_DD_fmt[i].fmt;
            getter = write_DD_fmt[i].getter;
            datatype = write_DD_fmt[i].datatype;
            datum_size = write_DD_fmt[i].datum_size;
            word_size = write_DD_fmt[i].word_size;
            datacount = 1;
            typesize = datum_size;
            break;
        }
    }
    if (getter == 0)
        return luaL_error(L, "qcd.ddpairs.write(): unknown precision");
    f_xml = QIO_string_create();
    QIO_string_set(f_xml, file_xml);
    s_xml = QIO_string_create();
    QIO_string_set(s_xml, source_xml);
    p_xml = QIO_string_create();
    QIO_string_set(p_xml, prop_xml);
    writer = qlua_qio_std_writer(fname, f_xml);
    if (writer == 0)
        return luaL_error(L, "qcd.ddpairs.write(): open error");
    dataP.P = QDP_expose_P(P->ptr);
    dataP.L = L;
    dataS.P = QDP_expose_P(S->ptr);
    dataS.L = L;
    for (js = 0; js < QDP_Ns; js++) {
        for (jc = 0; jc < QDP_Nc; jc++) {
            QIO_RecordInfo *rec_info;
            int rec_cl; /* This is how qio.h defines it */
            int *rec_lo, *rec_hi;

            dataS.js = js;
            dataS.jc = jc;
            if ((time_slice >= 0) && (time_slice < qDim[qRank - 1])) {
                rec_cl = QIO_HYPER;
                rec_lo = lo;
                rec_hi = hi;
            } else {
                rec_cl = QIO_FIELD;
                rec_lo = NULL;
                rec_hi = NULL;
            }
            rec_info = QIO_create_record_info(rec_cl, rec_lo, rec_hi, qRank,
                                              datatype, prec, QDP_Nc, QDP_Ns,
                                              typesize, datacount);
            if (QIO_write(writer, rec_info, s_xml, getter,
                          datum_size, word_size, &dataS))
                return luaL_error(L, "qcd.ddpairs.write(): source write error");
            QIO_destroy_record_info(rec_info);
            
            dataP.js = js;
            dataP.jc = jc;
            rec_info = QIO_create_record_info(QIO_FIELD, NULL, NULL, qRank,
                                              datatype, prec, QDP_Nc, QDP_Ns,
                                              typesize, datacount);
            if (QIO_write(writer, rec_info, p_xml, getter,
                          datum_size, word_size, &dataP))
                return luaL_error(L, "qcd.ddpairs.write(): source write prop");
            QIO_destroy_record_info(rec_info);
        }
    }

    QDP_reset_P(S->ptr);
    QDP_reset_P(P->ptr);
    QIO_string_destroy(p_xml);
    QIO_string_destroy(s_xml);
    QIO_string_destroy(f_xml);
    qlua_free(L, lo);
    qlua_free(L, hi);

    return 0;
}

static struct luaL_Reg fDDPAIRSio[] = {
    { "read",     ddpairs_read },
    { "write",    ddpairs_write },
    { NULL,       NULL }
};

int
init_ddpairs_io(lua_State *L)
{
    int i;

    lua_getglobal(L, qcdlib);
    lua_newtable(L);
    for (i = 0; fDDPAIRSio[i].name; i++) {
        lua_pushcfunction(L, fDDPAIRSio[i].func);
        lua_setfield(L, -2, fDDPAIRSio[i].name);
    }
    lua_setfield(L, -2, ddpairs_io);
    lua_pop(L, 1);
    return 0;
}

int
fini_ddpairs_io(lua_State *L)
{
    return 0;
}
