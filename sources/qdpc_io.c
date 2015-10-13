#include "modules.h"                                                 /* DEPS */
#include "qlua.h"                                                    /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "qdpc_io.h"                                                 /* DEPS */
#include "qio_utils.h"                                               /* DEPS */
#include "latcolmat.h"                                               /* DEPS */
#include "latcolvec.h"                                               /* DEPS */
#include "latcomplex.h"                                              /* DEPS */
#include "latdirferm.h"                                              /* DEPS */
#include "latdirprop.h"                                              /* DEPS */
#include "latint.h"                                                  /* DEPS */
#include "latrandom.h"                                               /* DEPS */
#include "latreal.h"                                                 /* DEPS */
#include "qvector.h"                                                 /* DEPS */
#ifdef HAS_GSL
#include "qmatrix.h"                                                 /* DEPS */
#endif

#include "qdp_df.h"
#if USE_Nc2
#include "qdp_f2.h"
#include "qdp_df2.h"
#endif
#if USE_Nc3
#include "qdp_f3.h"
#include "qdp_df3.h"
#endif
#if USE_NcN
#include "qdp_fn.h"
#include "qdp_dfn.h"
#endif
#include <string.h>

typedef struct {
    int master;
    int rank_stride;
} qdpc_ionode_ctrl_s;

static int qdpc_my_ionode_a(int node, void *arg) 
{
    qdpc_ionode_ctrl_s *ionode_ctrl = (qdpc_ionode_ctrl_s *)arg;
    int rank_stride = ionode_ctrl->rank_stride;
    return (node / rank_stride) * rank_stride;
}
static int qdpc_master_ionode_a(void *arg) 
{
    qdpc_ionode_ctrl_s *ionode_ctrl = (qdpc_ionode_ctrl_s *)arg;
    return ionode_ctrl->master;
}

typedef struct {
    QDP_Reader *ptr;
    qdpc_ionode_ctrl_s *ionode_ctrl;
} mReader;

typedef struct {
    QDP_Writer *ptr;
    qdpc_ionode_ctrl_s *ionode_ctrl;
} mWriter;

/* this magic used to identify global matrices in QIO */
static const char qio_matrix_magic[] = "480C94A9-F941-4FD5-9882-3E70B004AACF";

static const char qdp_io[] = "qdpc";
static const char ReaderName[] = "qcd.qdpc.reader";
static const char WriterName[] = "qcd.qdpc.writer";

static mReader *q_checkReader(lua_State *L, int idx, mLattice *S);
static mWriter *q_checkWriter(lua_State *L, int idx, mLattice *S);

/* check for a live QIO object */
static void check_reader(lua_State *L, mReader *b);
static void check_writer(lua_State *L, mWriter *b);

/* get QIO record precision */
static int
get_prec(QDP_Reader *qr)
{
  QIO_RecordInfo *ri = QIO_create_record_info(0, NULL, NULL, 0, "", "", 0, 0, 0, 0);
  QDP_String *md = QDP_string_create();
  int prec;

  QDP_read_qio_record_info(qr, ri, md);
  prec = *QIO_get_precision(ri);
  QIO_destroy_record_info(ri);
  QDP_string_destroy(md);

  return prec;
}

/* plain lattice types */
#define T_QTYPE       QDP_Int
#define QLUA_NAME(x)  x ## LatInt
#define X_ID(x)       x ## I
#undef X_DF
#include "qdpc_io-x.c"                                               /* DEPS */

#define T_QTYPE       QDP_RandomState
#define QLUA_NAME(x)  x ## LatRandom
#define X_ID(x)       x ## S
#undef X_DF
#include "qdpc_io-x.c"                                               /* DEPS */

#define T_QTYPE       QDP_D_Real
#define QLUA_NAME(x)  x ## LatReal
#define X_ID(x)       x ## R
#define X_ID2(a,b)    a ## R ## b
#define X_ID3(a,b)    a ## R ## b ## R
#define T_sTYPE       QDP_F_Real
#define SopL(n)       QDP_F_ ## n ## _R_L
#define Sop(n)        QDP_F_ ## n ## _R
#define FtoD(d,s,a)   QDP_DF_R_eq_R((d),(s),(a))
#define X_DF
#include "qdpc_io-x.c"                                               /* DEPS */

#define T_QTYPE       QDP_D_Complex
#define QLUA_NAME(x)  x ## LatComplex
#define X_ID(x)       x ## C
#define X_ID2(a,b)    a ## C ## b
#define X_ID3(a,b)    a ## C ## b ## C
#define T_sTYPE       QDP_F_Complex
#define SopL(n)       QDP_F_ ## n ## _C_L
#define Sop(n)        QDP_F_ ## n ## _C
#define FtoD(d,s,a)   QDP_DF_C_eq_C((d),(s),(a))
#define X_DF
#include "qdpc_io-x.c"                                               /* DEPS */

#if USE_Nc2 || USE_Nc3 || USE_NcN
/* colored lattice types */
#define QT(a)         a ## ColorVector
#define QTx(a,b)      a ## ColorVector ## b
#define QN(a,b)       a ## LatColVec ## b
#define QA(p)         p ## V
#define QAx(a,b)      a ## V ## b
#define QAy(a,b)      a ## V ## b ## V
#include "qdpc_io-y.c"                                               /* DEPS */

#define QT(a)         a ## ColorMatrix
#define QTx(a,b)      a ## ColorMatrix ## b
#define QN(a,b)       a ## LatColMat ## b
#define QA(p)         p ## M
#define QAx(a,b)      a ## M ## b
#define QAy(a,b)      a ## M ## b ## M
#include "qdpc_io-y.c"                                               /* DEPS */

#define QT(a)         a ## DiracFermion
#define QTx(a,b)      a ## DiracFermion ## b
#define QN(a,b)       a ## LatDirFerm ## b
#define QA(p)         p ## D
#define QAx(a,b)      a ## D ## b
#define QAy(a,b)      a ## D ## b ## D
#include "qdpc_io-y.c"                                               /* DEPS */

#define QT(a)         a ## DiracPropagator
#define QTx(a,b)      a ## DiracPropagator ## b
#define QN(a,b)       a ## LatDirProp ## b
#define QA(p)         p ## P
#define QAx(a,b)      a ## P ## b
#define QAy(a,b)      a ## P ## b ## P
#include "qdpc_io-y.c"                                               /* DEPS */
#endif /* use any colors */

/* QDPC reader */
static int
qdpc_r_fmt(lua_State *L)
{
    mReader *b = q_checkReader(L, 1, NULL);
    char fmt[72];
    
    if ( b->ptr )
        sprintf(fmt, "qdpc.Reader(%p)", b->ptr);
    else
        sprintf(fmt, "qdpc.Reader(closed)");

    lua_pushstring(L, fmt);

    return 1;
}

static int
qdpc_r_gc(lua_State *L)
{
    mReader *b = q_checkReader(L, 1, NULL);

    if (b->ptr)
        QDP_close_read(b->ptr);
    b->ptr = 0;
    if (NULL != b->ionode_ctrl) 
        qlua_free(L, b->ionode_ctrl);
    b->ionode_ctrl = NULL;
    
    return 0;
}

static void
check_reader(lua_State *L, mReader *b)
{
    if (b->ptr == 0)
        luaL_error(L, "closed QDP/C reader");
}

static int
qdpc_r_skip(lua_State *L)
{
    mReader *reader = q_checkReader(L, 1, NULL);
    int status;
    int rcount;

    check_reader(L, reader);

    CALL_QDP(L);
    status = QDP_next_record(reader->ptr);
    if (status == 0) {
        lua_pushboolean(L, 1);
        rcount = 1;
    } else {
        rcount = 0;
    }

    return rcount;
}

static int
qdpc_r_close(lua_State *L)
{
    mReader *b = q_checkReader(L, 1, NULL);
    
    check_reader(L, b);

    QDP_close_read(b->ptr);
    b->ptr = 0;
    
    return 0;
}

static int
qdpc_r_iv(lua_State *L)
{
    mReader *reader = q_checkReader(L, 1, NULL);
    int len = luaL_checkint(L, 2);
    mVecInt *v = qlua_newVecInt(L, len);
    QDP_String *info;
    
    check_reader(L, reader);

    CALL_QDP(L);

    info = QDP_string_create();
    if (len <= 0)
        goto error;

    if (QDP_vread_i(reader->ptr, info, v->val, v->size) != 0)
        goto error;

    lua_pushstring(L, QDP_string_ptr(info));
    QDP_string_destroy(info);
    return 2;
error:
    QDP_string_destroy(info);
    return luaL_error(L, "qdpc read error");
}

static int
qdpc_r_rv(lua_State *L)
{
    mReader *reader = q_checkReader(L, 1, NULL);
    int len = luaL_checkint(L, 2);
    mVecReal *v = qlua_newVecReal(L, len);
    QDP_String *info;
    
    check_reader(L, reader);

    CALL_QDP(L);

    info = QDP_string_create();
    if (len <= 0)
        goto error;

    if (QDP_D_vread_r(reader->ptr, info, v->val, v->size) != 0)
        goto error;

    lua_pushstring(L, QDP_string_ptr(info));
    QDP_string_destroy(info);
    return 2;
error:
    QDP_string_destroy(info);
    return luaL_error(L, "qdpc read error");
}

static int
qdpc_r_cv(lua_State *L)
{
    mReader *reader = q_checkReader(L, 1, NULL);
    int len = luaL_checkint(L, 2);
    mVecComplex *v = qlua_newVecComplex(L, len);
    QDP_String *info;
    
    check_reader(L, reader);

    CALL_QDP(L);

    info = QDP_string_create();
    if (len <= 0)
        goto error;

    if (QDP_D_vread_c(reader->ptr, info, v->val, v->size) != 0)
        goto error;

    lua_pushstring(L, QDP_string_ptr(info));
    QDP_string_destroy(info);
    return 2;
error:
    QDP_string_destroy(info);
    return luaL_error(L, "qdpc read error");
}

typedef struct {
    QLA_D_Real **qla_arr;
    int len;
} putR_func_arg_s; 
/* qdpc_r_genericR_put_F|D 
   * read a sequence of float|double
   * ignore "count" 
*/
static void 
qdpc_r_genericR_put_F(char *buf, size_t index, int count, void *arg_)
{
    putR_func_arg_s *arg = arg_;
//    if (arg->len != count)
//        printf("arg->len =%d != %d=count\n", arg->len, count);
    for (int i = 0 ; i < arg->len ; i++)
        arg->qla_arr[i][index] = ((float *)buf)[i];
}
static void 
qdpc_r_genericR_put_D(char *buf, size_t index, int count, void *arg_)
{
    putR_func_arg_s *arg = arg_;
//    if (arg->len != count)
//        printf("arg->len =%d != %d=count\n", arg->len, count);
    for (int i = 0 ; i < arg->len ; i++)
        arg->qla_arr[i][index] = ((double *)buf)[i];
}

static int
qdpc_r_GenericReal(lua_State *L) 
{
    /* read any lime record interpreting it as a sequence of LatticeReal (determi) 
       return a table (list) of LatticeReal
       determine the length of 
     */
    int status = 0;
    const char *errstr = NULL;

    mReader *reader = q_checkReader(L, 1, NULL);
    int Sidx;


    QDP_D_Real **qdp_arr = NULL;
    QLA_D_Real **qla_arr = NULL;
    QIO_String *rec_xml  = NULL;
    char prec = '\0' ;   /* no default value */
    int datacount   = -1,
        typesize    = -1,
        wordsize    = -1;

    void (*put_f)(char *, size_t, int, void *);
    putR_func_arg_s put_arg;

    if (qlua_checkopt_table(L, 2)) {
        if (qlua_tabpushopt_key(L, 2, "precision")) {
            prec = qlua_qio_file_precision(L, -1);
            lua_pop(L, 1);
        }
    }

    check_reader(L, reader);
    qlua_ObjLattice(L, 1);
    Sidx = lua_gettop(L);
    
    CALL_QDP(L);

    QIO_Reader *qio_r = QDP_reader_get_qio(reader->ptr);
    QIO_RecordInfo rec_info;
    if (NULL == (rec_xml = QIO_string_create())) {
        errstr = "memory error";
        goto clearerr_1;
    }
    if (0 != (status =
            QIO_read_record_info(qio_r, &rec_info, rec_xml))) {
        errstr = "cannot read record info";
        goto clearerr_1;
    }
    if ('F' != prec && 'D' != prec) {
        if (!QIO_defined_precision(&rec_info)) {
            luaL_error(L, "precision not defined");
        }
        prec = *QIO_get_precision(&rec_info);
    }
    switch (prec) {
    case 'F' :  {
        wordsize = 4;       
        put_f    = qdpc_r_genericR_put_F;
        break; }
    case 'D' :  {
        wordsize = 8;       
        put_f    = qdpc_r_genericR_put_D;
        break;
    }
    default:    return luaL_error(L, "unknown precision");
    }
    if (!QIO_defined_datacount(&rec_info) || !QIO_defined_typesize(&rec_info)) {
        errstr = "datacount or typesize not defined";
        goto clearerr_1;
    }
    datacount   = QIO_get_datacount(&rec_info);
    typesize    = QIO_get_typesize(&rec_info);
    if (0 != (datacount * typesize) % wordsize) {
        errstr = "incompatible datacount, typesize and precision";
        goto clearerr_1;
    }
    int lenR        = (datacount * typesize) / wordsize;
//    printf("prec='%c'   datacount=%d   typesize=%d   wordsize=%d\n", prec, datacount, typesize, wordsize);

    /* init fields, prepare for import */
    qdp_arr = qlua_malloc(L, lenR * sizeof(QDP_D_Real *));
    qla_arr = qlua_malloc(L, lenR * sizeof(QLA_D_Real *));
    if (NULL == qdp_arr || NULL == qla_arr) {
        errstr = "memory_error";
        goto clearerr_2;
    }

    lua_createtable(L, lenR, 0);
    for (int i = 0 ; i < lenR ; i++) {
        mLatReal *ui = qlua_newZeroLatReal(L, Sidx);
        if (NULL == ui) {
            errstr = "memory_error";
            goto clearerr_2;
        }
        qdp_arr[i] = ui->ptr;
        qla_arr[i] = QDP_D_expose_R(qdp_arr[i]);
        lua_rawseti(L, -2, i + 1); /* pops the value */
    }

    /* import qio data */
    put_arg.len     = lenR;
    put_arg.qla_arr = qla_arr;
    if (0 != (status = QIO_read(qio_r, &rec_info, rec_xml, put_f, 
                    datacount * typesize, wordsize, &put_arg))) {
        errstr = "qio read failed";
        goto clearerr_2;
    }
    
    lua_pushstring(L, QIO_string_ptr(rec_xml));

    /* cleanup */
    for (int i = 0 ; i < lenR ; i++)
        QDP_D_reset_R(qdp_arr[i]);
    qlua_free(L, qdp_arr);
    qlua_free(L, qla_arr);
    QIO_string_destroy(rec_xml);

    return 2;       /* { (table of)field(s), rec_xml } */

    /* error cleanup */
clearerr_3:
    for (int i = 0 ; i < lenR ; i++) {
        QDP_D_reset_R(qdp_arr[i]);
        QDP_D_destroy_R(qdp_arr[i]);
    }
clearerr_2:
    if (NULL != qla_arr) qlua_free(L, qla_arr);
    if (NULL != qdp_arr) qlua_free(L, qdp_arr);
clearerr_1:
    if (NULL != rec_xml) QIO_string_destroy(rec_xml);

    return luaL_error(L, errstr);
}



#ifdef HAS_GSL
static int
qdpc_r_rm(lua_State *L)
{
    mReader *reader = q_checkReader(L, 1, NULL);
    int l_size = luaL_checkint(L, 2);
    int r_size = luaL_checkint(L, 3);
    mMatReal *v = qlua_newMatReal(L, l_size, r_size);
    QLA_Int f_size[2];
    QDP_String *info;
    int i, j;

    check_reader(L, reader);

    CALL_QDP(L);

    info = QDP_string_create();
    if ((l_size <= 0) || (r_size <= 0))
        goto error;

    if (QDP_vread_i(reader->ptr, info, f_size, 2) != 0)
        goto error;
    if ((f_size[0] != l_size) ||
        (f_size[1] != r_size))
        goto error;
    if (strcmp(QDP_string_ptr(info), qio_matrix_magic) != 0)
        goto error;
    {
        QLA_D_Real d[l_size * r_size];
        if (QDP_D_vread_r(reader->ptr, info, d, l_size * r_size) != 0) {
            goto error;
        }
        /* must only agree with the writer below */
        for (j = 0; j < r_size; j++) {
            for (i = 0; i < l_size; i++) {
                gsl_matrix_set(v->m, i, j, d[i + l_size * j]);
            }
        }
        lua_pushstring(L, QDP_string_ptr(info));
        QDP_string_destroy(info);
        return 2;
    }
error:
    QDP_string_destroy(info);
    return luaL_error(L, "qdpc read error");
}

static int
qdpc_r_cm(lua_State *L)
{
    mReader *reader = q_checkReader(L, 1, NULL);
    int l_size = luaL_checkint(L, 2);
    int r_size = luaL_checkint(L, 3);
    mMatComplex *v = qlua_newMatComplex(L, l_size, r_size);
    QLA_Int f_size[2];
    QDP_String *info;
    int i, j;

    check_reader(L, reader);

    CALL_QDP(L);

    info = QDP_string_create();
    if ((l_size <= 0) || (r_size <= 0))
        goto error;

    if (QDP_vread_i(reader->ptr, info, f_size, 2) != 0)
        goto error;
    if ((f_size[0] != l_size) ||
        (f_size[1] != r_size))
        goto error;
    if (strcmp(QDP_string_ptr(info), qio_matrix_magic) != 0)
        goto error;
    {
        QLA_D_Complex d[l_size * r_size];
        if (QDP_D_vread_c(reader->ptr, info, d, l_size * r_size) != 0) {
            goto error;
        }
        /* must only agree with the writer below */
        for (j = 0; j < r_size; j++) {
            for (i = 0; i < l_size; i++) {
                QLA_D_Complex z = d[i + l_size * j];
                gsl_complex zz;
                
                GSL_SET_COMPLEX(&zz, QLA_real(z), QLA_imag(z));
                gsl_matrix_complex_set(v->m, i, j, zz);
            }
        }
        lua_pushstring(L, QDP_string_ptr(info));
        QDP_string_destroy(info);
        return 2;
    }
error:
    QDP_string_destroy(info);
    return luaL_error(L, "qdpc read error");
}
#endif

static const struct luaL_Reg mtReader[] = {
    { "__tostring",       qdpc_r_fmt               },
    { "__gc",             qdpc_r_gc                },
    { "__newindex",       qlua_nowrite             },
    { "skip",             qdpc_r_skip              },
    { "close",            qdpc_r_close             },
    { "Int",              qdpc_r_I                 },
    { "RandomState",      qdpc_r_S                 },
    { "Real",             qdpc_r_R                 },
    { "Complex",          qdpc_r_C                 },
#if USE_Nc2 || USE_Nc3 || USE_NcN
    { "ColorVector",      qdpc_r_ColorVector       },
    { "generic_Real",     qdpc_r_GenericReal       },
    { "ColorVectorN",     qdpc_r_ColorVectorN      },
    { "ColorMatrix",      qdpc_r_ColorMatrix       },
    { "ColorMatrixN",     qdpc_r_ColorMatrixN      },
    { "DiracFermion",     qdpc_r_DiracFermion      },
    { "DiracFermionN",    qdpc_r_DiracFermionN     },
    { "DiracPropagator",  qdpc_r_DiracPropagator   },
    { "DiracPropagatorN", qdpc_r_DiracPropagatorN  },
#endif /* use any colors */
    { "int_vector",       qdpc_r_iv                },
    { "real_vector",      qdpc_r_rv                },
    { "complex_vector",   qdpc_r_cv                },
#ifdef HAS_GSL
    { "real_matrix",      qdpc_r_rm                },
    { "complex_matrix",   qdpc_r_cm                },
#endif
    /* "lattice" */
    /* "a-type" */
    { NULL,               NULL                     }
};

static mReader *
q_checkReader(lua_State *L, int idx, mLattice *S)
{
    void *v = qlua_checkLatticeType(L, idx, qReader, ReaderName);
    
    if (S) {
        mLattice *S1 = qlua_ObjLattice(L, idx);
        if (S1->id != S->id)
            luaL_error(L, "%s on a wrong lattice", ReaderName);
        lua_pop(L, 1);
    }

    return (mReader *)v;
}

static mReader *
q_newReader(lua_State *L, QDP_Reader *reader, qdpc_ionode_ctrl_s *ionode_ctrl)
{
    mReader *h = lua_newuserdata(L, sizeof (mReader));
    h->ptr          = reader;
    h->ionode_ctrl  = ionode_ctrl;

    qlua_createLatticeTable(L, 1, mtReader, qReader, ReaderName);
    lua_setmetatable(L, -2);

    return h;
}

/* QDP/C reader constructors
 *
 * qcd.qdpc.reader(L, file_name)
 */
static int
q_qdpc_reader (lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);
    const char *name = luaL_checkstring(L, 2);
    int has_fs = 0,
        ionode_rank_stride = 1,
        ionode_master = 0;
    qdpc_ionode_ctrl_s *ionode_ctrl = NULL;
    QIO_Filesystem fs;
    if (qlua_checkopt_table(L, 3)) {
        if (qlua_tabpushopt_key(L, 3, "rank_stride")) {
            ionode_rank_stride = luaL_checkint(L, -1);
            has_fs  = 1;
            lua_pop(L, 1);
        }
        if (qlua_tabpushopt_key(L, 3, "master")) {
            ionode_master = luaL_checkint(L, -1);
            has_fs  = 1;
            lua_pop(L, 1);
        }
    }
    
    QDP_String *xml = 0;
    QDP_Reader *reader = 0;
    int rcount;

    /* first collect garbage */
    CALL_QDP(L);

    /* go through the motions of opening a QDP reader */
    xml = QDP_string_create();
    if (has_fs) {
        fs.my_io_node_a         = qdpc_my_ionode_a;
        fs.master_io_node_a     = qdpc_master_ionode_a;
        ionode_ctrl             = qlua_malloc(L, sizeof(qdpc_ionode_ctrl_s));
        if (NULL == ionode_ctrl) 
            luaL_error(L, "memory error");
        ionode_ctrl->master     = ionode_master;
        ionode_ctrl->rank_stride= ionode_rank_stride;
        fs.arg                  = ionode_ctrl;
        /* TODO perhaps set only those overridden by lua params */
        reader  = QDP_open_read_general_L(S->lat, xml, (char *)name, &fs, NULL);
    } else {
        reader  = QDP_open_read_L(S->lat, xml, (char *)name); /* [sic] */
    }

    /* convert QDP results to LUA */
    if (reader == 0) {
        /* reader is not opened - return no results */
        rcount = 0;
    } else {
        /* we have a reader - return (mReader, file_xml) */
        q_newReader(L, reader, ionode_ctrl);
        lua_pushstring(L, QDP_string_ptr(xml));
        rcount = 2;
    }
    /* free QDP_string now */
    QDP_string_destroy(xml);

    return rcount;
}

/* QDPC writer */

static void
check_writer(lua_State *L, mWriter *b)
{
    if (b->ptr == 0)
        luaL_error(L, "closed QDP/C writer");
}

static int
qdpc_w_fmt(lua_State *L)
{
    mWriter *b = q_checkWriter(L, 1, NULL);
    char fmt[72];
    
    if (b->ptr)
        sprintf(fmt, "qdpc.Writer(%p)", b->ptr);
    else
        sprintf(fmt, "qdpc.Writer(closed)");

    lua_pushstring(L, fmt);

    return 1;
}

static int
qdpc_w_gc(lua_State *L)
{
    mWriter *b = q_checkWriter(L, 1, NULL);

    if (b->ptr)
        QDP_close_write(b->ptr);
    b->ptr = 0;
    if (NULL != b->ionode_ctrl) 
        qlua_free(L, b->ionode_ctrl);
    b->ionode_ctrl = NULL;
    
    return 0;
}

static int
qdpc_w_close(lua_State *L)
{
    mWriter *b = q_checkWriter(L, 1, NULL);

    check_writer(L, b);
    
    QDP_close_write(b->ptr);
    b->ptr = 0;

    return 0;
}

static int
qdpc_w_iv(lua_State *L)
{
    mWriter *writer = q_checkWriter(L, 1, NULL);
    mVecInt *v = qlua_checkVecInt(L, 2);
    const char *info = luaL_checkstring(L, 3);
    QDP_String *xml;
    
    check_writer(L, writer);

    CALL_QDP(L);

    xml = QDP_string_create();
    QDP_string_set(xml, (char *)info); /* [ sic ] */

    if (QDP_vwrite_i(writer->ptr, xml, v->val, v->size) == 0) {
        lua_pushboolean(L, 1);
        QDP_string_destroy(xml);
        return 1;
    } else {
        QDP_string_destroy(xml);
        return luaL_error(L, "qdpc write error");
    }
}

static int
qdpc_w_rv(lua_State *L)
{
    mWriter *writer = q_checkWriter(L, 1, NULL);
    mVecReal *v = qlua_checkVecReal(L, 2);
    const char *info = luaL_checkstring(L, 3);
    QDP_String *xml;
    
    check_writer(L, writer);

    CALL_QDP(L);

    xml = QDP_string_create();
    QDP_string_set(xml, (char *)info); /* [ sic ] */

    if (QDP_D_vwrite_r(writer->ptr, xml, v->val, v->size) == 0) {
        lua_pushboolean(L, 1);
        QDP_string_destroy(xml);
        return 1;
    } else {
        QDP_string_destroy(xml);
        return luaL_error(L, "qdpc write error");
    }
}

static int
qdpc_w_cv(lua_State *L)
{
    mWriter *writer = q_checkWriter(L, 1, NULL);
    mVecComplex *v = qlua_checkVecComplex(L, 2);
    const char *info = luaL_checkstring(L, 3);
    QDP_String *xml;
    
    check_writer(L, writer);

    CALL_QDP(L);

    xml = QDP_string_create();
    QDP_string_set(xml, (char *)info); /* [ sic ] */

    if (QDP_D_vwrite_c(writer->ptr, xml, v->val, v->size) == 0) {
        lua_pushboolean(L, 1);
        QDP_string_destroy(xml);
        return 1;
    } else {
        QDP_string_destroy(xml);
        return luaL_error(L, "qdpc write error");
    }
}

#ifdef HAS_GSL
static int
qdpc_w_rm(lua_State *L)
{
    mWriter *writer = q_checkWriter(L, 1, NULL);
    mMatReal *v = qlua_checkMatReal(L, 2);
    const char *info = luaL_checkstring(L, 3);
    int l_size = v->l_size;
    int r_size = v->r_size;
    QLA_Int f_size[2];
    QDP_String *xml;
    int i, j;

    check_writer(L, writer);

    f_size[0] = l_size;
    f_size[1] = r_size;
    QLA_D_Real d[l_size * r_size];
    /* must only agree with the reader above */
    for (j = 0; j < r_size; j++) {
        for (i = 0; i < l_size; i++) {
            d[i + l_size * j] = gsl_matrix_get(v->m, i, j);
        }
    }

    CALL_QDP(L);

    xml = QDP_string_create();
    QDP_string_set(xml, (char *)qio_matrix_magic); /* [ sic ] */
    if (QDP_vwrite_i(writer->ptr, xml, f_size, 2) != 0)
        goto error;
    QDP_string_set(xml, (char *)info); /* [ sic ] */
    if (QDP_D_vwrite_r(writer->ptr, xml, d, l_size * r_size) != 0)
        goto error;
    QDP_string_destroy(xml);
    lua_pushboolean(L, 1);
    return 1;
error:
    QDP_string_destroy(xml);
    return luaL_error(L, "qdpc write error");
}

static int
qdpc_w_cm(lua_State *L)
{
    mWriter *writer = q_checkWriter(L, 1, NULL);
    mMatComplex *v = qlua_checkMatComplex(L, 2);
    const char *info = luaL_checkstring(L, 3);
    int l_size = v->l_size;
    int r_size = v->r_size;
    QLA_Int f_size[2];
    QDP_String *xml;
    int i, j;

    check_writer(L, writer);

    f_size[0] = l_size;
    f_size[1] = r_size;
    QLA_D_Complex *d = qlua_malloc(L, l_size * r_size * sizeof (QLA_D_Complex));
    /* must only agree with the reader above */
    for (j = 0; j < r_size; j++) {
        for (i = 0; i < l_size; i++) {
            gsl_complex z = gsl_matrix_complex_get(v->m, i, j);
            QLA_real(d[i + l_size * j]) = GSL_REAL(z);
            QLA_imag(d[i + l_size * j]) = GSL_IMAG(z);
        }
    }

    CALL_QDP(L);

    xml = QDP_string_create();
    QDP_string_set(xml, (char *)qio_matrix_magic); /* [ sic ] */
    if (QDP_vwrite_i(writer->ptr, xml, f_size, 2) != 0)
        goto error;
    QDP_string_set(xml, (char *)info); /* [ sic ] */
    if (QDP_D_vwrite_c(writer->ptr, xml, d, l_size * r_size) != 0)
        goto error;
    QDP_string_destroy(xml);
    lua_pushboolean(L, 1);
        qlua_free(L, d);
    return 1;
error:
        qlua_free(L, d);
    QDP_string_destroy(xml);
    return luaL_error(L, "qdpc write error");
}
#endif

static const struct luaL_Reg mtWriter[] = {
    { "__tostring",       qdpc_w_fmt             },
    { "__gc",             qdpc_w_gc              },
    { "__newindex",       qlua_nowrite           },
    { "close",            qdpc_w_close           },
    { "Int",              qdpc_w_I               },
    { "RandomState",      qdpc_w_S               },
    { "Real",             qdpc_w_R               },
    { "Complex",          qdpc_w_C               },
#if USE_Nc2 || USE_Nc3 || USE_NcN
    { "ColorVector",      qdpc_w_ColorVector     },
    { "ColorMatrix",      qdpc_w_ColorMatrix     },
    { "DiracFermion",     qdpc_w_DiracFermion    },
    { "DiracPropagator",  qdpc_w_DiracPropagator },
#endif /* use any colors */
    { "int_vector",       qdpc_w_iv              },
    { "real_vector",      qdpc_w_rv              },
    { "complex_vector",   qdpc_w_cv              },
#ifdef HAS_GSL
    { "real_matrix",      qdpc_w_rm              },
    { "complex_matrix",   qdpc_w_cm              },
#endif
    /* "lattice" */
    /* "a-type" */
    { NULL,               NULL                   }
};

static mWriter *
q_checkWriter(lua_State *L, int idx, mLattice *S)
{
    void *v = qlua_checkLatticeType(L, idx, qWriter, WriterName);
    
    if (S) {
        mLattice *S1 = qlua_ObjLattice(L, idx);
        if (S1->id != S->id)
            luaL_error(L, "%s on a wrong lattice", WriterName);
        lua_pop(L, 1);
    }

    return (mWriter *)v;
}

static mWriter *
q_newWriter(lua_State *L, QDP_Writer *writer, qdpc_ionode_ctrl_s *ionode_ctrl)
{
    mWriter *h = lua_newuserdata(L, sizeof (mWriter));
    h->ptr          = writer;
    h->ionode_ctrl  = ionode_ctrl;

    qlua_createLatticeTable(L, 1, mtWriter, qWriter, WriterName);
    lua_setmetatable(L, -2);

    return h;
}

/* QDP/C writer constructors
 *
 *  qcd.qdpc.writer(L, file_name, file_info)
 *  qcd.qdpc.writer(L, file_name, file_info, "single")
 *  qcd.qdpc.writer(L, file_name, file_info, "multi")
 */
static int
q_qdpc_writer(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);
    const char *name = luaL_checkstring(L, 2);
    const char *info = luaL_checkstring(L, 3);
    int volfmt = qlua_qio_volume_format(L, 4, lua_gettop(L));
    QDP_String *xml = 0;
    QDP_Writer *writer = 0;
    int rcount;

    int has_fs = 0,
        ionode_rank_stride = INT_MAX,
        ionode_master = 0;
    qdpc_ionode_ctrl_s *ionode_ctrl = NULL;
    QIO_Filesystem fs;
    if (qlua_checkopt_table(L, 5)) {
        if (qlua_tabpushopt_key(L, 5, "rank_stride")) {
            ionode_rank_stride  = luaL_checkint(L, -1);
            has_fs              = 1;
            lua_pop(L, 1);
        }
        if (qlua_tabpushopt_key(L, 5, "master")) {
            ionode_master   = luaL_checkint(L, -1);
            has_fs          = 1;
            lua_pop(L, 1);
        }
    }

    /* first collect garbage */
    CALL_QDP(L);
    
    /* open the QDP writer */
    xml = QDP_string_create();
    QDP_string_set(xml, (char *)info); /* [ sic ] */
    if (has_fs) {
        fs.my_io_node_a         = qdpc_my_ionode_a;
        fs.master_io_node_a     = qdpc_master_ionode_a;
        ionode_ctrl             = qlua_malloc(L, sizeof(qdpc_ionode_ctrl_s));
        if (NULL == ionode_ctrl) 
            luaL_error(L, "memory error");
        ionode_ctrl->master     = ionode_master;
        ionode_ctrl->rank_stride= ionode_rank_stride;
        fs.arg                  = ionode_ctrl;
        /* TODO perhaps set only those overridden by lua params */
        writer = QDP_open_write_general_L(S->lat, xml, (char *)name, volfmt, &fs, NULL); /* [ sic ] */
    } else {
        writer = QDP_open_write_L(S->lat, xml, (char *)name, volfmt); /* [ sic ] */
    }

    /* convert QDP results into LUA values */
    if (writer == 0) {
        /* failure - zero results */
        rcount = 0;
    } else {
        /* success -- return writer */
        q_newWriter(L, writer, ionode_ctrl);
        rcount = 1;
    }
    QDP_string_destroy(xml);

    return rcount;
}

/* names and routines for qcd.qdpc table */
static const struct luaL_Reg fQDPio[] = {
    { "Reader",   q_qdpc_reader},
    { "Writer",   q_qdpc_writer},
    { NULL,       NULL }
};

int
init_qdpc_io(lua_State *L)
{
    lua_getglobal(L, qcdlib);
    lua_newtable(L);
    luaL_register(L, NULL, fQDPio);
    lua_setfield(L, -2, qdp_io);
    lua_pop(L, 1);
    
    return 0;
}

void
fini_qdpc_io(void)
{
}
