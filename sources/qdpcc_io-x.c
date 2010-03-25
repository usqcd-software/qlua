
typedef struct {
    int         nc;
    lua_State  *L;
    void       *P; /* QLA_DiracPropagator *P */
} Qs(QDPCCArgs);

static void
Qs(float_P_put)(char *buf, size_t index, int count, void *v_arg)
{
    Qs(QDPCCArgs) *arg = v_arg;
#if QNc == 'N'
    typedef QLA_DN_DiracPropagator(arg->nc, Ptype);
#else
    typedef Qx(QLA_D,_DiracPropagator) Ptype;
#endif
    Ptype *P = arg->P;
    Ptype *dst = &P[index];
    QLA_F_Complex *src = (void *)buf;
    int ic, jc, is, js;

    if (count != 1)
        luaL_error(arg->L, "qcd.qdpcc.read_prop() count != 1");

    for (is = 0; is < QDP_Ns; is++) {
        for (js = 0; js < QDP_Ns; js++) {
            for (ic = 0; ic < arg->nc; ic++) {
                for (jc = 0; jc < arg->nc; jc++, src++) {
                    QLA_c_eq_r_plus_ir(QLA_elem_P(*dst, ic, is, jc, js),
                                       QLA_real(*src), QLA_imag(*src));
                               
                }
            }
        }
    }
}

static void
Qs(double_P_put)(char *buf, size_t index, int count, void *v_arg)
{
    Qs(QDPCCArgs) *arg = v_arg;
#if QNc == 'N'
    typedef QLA_DN_DiracPropagator(arg->nc, Ptype);
#else
    typedef Qx(QLA_D,_DiracPropagator) Ptype;
#endif
    Ptype *P = arg->P;
    Ptype *dst = &P[index];
    QLA_D_Complex *src = (void *)buf;
    int ic, jc, is, js;

    if (count != 1)
        luaL_error(arg->L, "qcd.qdpcc.read_prop() count != 1");

    for (is = 0; is < QDP_Ns; is++) {
        for (js = 0; js < QDP_Ns; js++) {
            for (ic = 0; ic < arg->nc; ic++) {
                for (jc = 0; jc < arg->nc; jc++, src++) {
                    QLA_c_eq_r_plus_ir(QLA_elem_P(*dst, ic, is, jc, js),
                                       QLA_real(*src), QLA_imag(*src));
                               
                }
            }
        }
    }
}

static int
Qs(cc_read_P_)(lua_State *L, int off, int nc)
{
#if QNc == 'N'
    typedef QLA_DN_DiracPropagator(nc, Ptype);
    typedef QLA_DN_DiracPropagator(nc, Dtype);
    typedef QLA_FN_DiracPropagator(nc, Ftype);
#else
    typedef Qx(QLA_D,_DiracPropagator) Ptype;
    typedef Qx(QLA_D,_DiracPropagator) Dtype;
    typedef Qx(QLA_F,_DiracPropagator) Ftype;
#endif
    typedef void Putter(char *buf, size_t index, int count, void *arg);

    mLattice         *S = qlua_checkLattice(L, 1 + off);
    const char       *fname = luaL_checkstring(L, 2 + off);
    Qs(mLatDirProp)  *prop = Qs(qlua_newLatDirProp)(L, 1 + off, nc);
    QIO_Reader       *reader = 0;
    QIO_String       *file_xml = 0;
    QIO_String       *record_xml = 0;
    QIO_RecordInfo   *record_info = 0;
    Putter           *putter = 0;
    size_t            datum_size = 0;
    int               word_size = 0;
    Qs(QDPCCArgs)     data;
    int               prec;

    CALL_QDP(L);
    file_xml = QIO_string_create();
    record_xml = QIO_string_create();
    reader = qlua_qio_std_reader(S, fname, file_xml);
    if (reader == 0)
        return luaL_error(L, "qcd.qdpcc.read_prop() open error");
    record_info = QIO_create_record_info(0, NULL, NULL, 0, "", "", 0, 0, 0, 0);
    if (QIO_read_record_info(reader, record_info, record_xml))
        return luaL_error(L, "qcd.qdpcc.read_prop() read info error");
    prec = *QIO_get_precision(record_info);
    switch (prec) {
    case 'F':
        putter = Qs(float_P_put);
        word_size = sizeof (QLA_F_Real);
        datum_size = sizeof (Ftype);
        break;
    case 'D':
        putter = Qs(double_P_put);
        word_size = sizeof (QLA_D_Real);
        datum_size = sizeof (Dtype);
        break;
    default:
        return luaL_error(L, "unsupported file precision");
    }
    data.nc = nc;
    data.P = Qx(QDP_D,_expose_P)(prop->ptr);
    data.L = L;
    if (QIO_read(reader, record_info, record_xml,
                 putter, datum_size, word_size, &data))
        luaL_error(L, "qcd.qdpcc.read_prop() read data error");
    Qx(QDP_D,_reset_P)(prop->ptr);
    QIO_destroy_record_info(record_info);
    QIO_close_read(reader);
    lua_pushstring(L, QIO_string_ptr(file_xml));
    lua_pushstring(L, QIO_string_ptr(record_xml));
    QIO_string_destroy(file_xml);
    QIO_string_destroy(record_xml);

    return 3;
}

static void
Qs(float_P_get)(char *buf, size_t index, int count, void *v_arg)
{
    Qs(QDPCCArgs) *arg = v_arg;
#if QNc == 'N'
    typedef QLA_DN_DiracPropagator(arg->nc, Ptype);
#else
    typedef Qx(QLA_D,_DiracPropagator) Ptype;
#endif
    Ptype *P = arg->P;
    Ptype *src = &P[index];
    QLA_F_Complex *dst = (void *)buf;
    int ic, jc, is, js;

    if (count != 1)
        luaL_error(arg->L, "qcd.qdpcc.write_prop() count != 1");

    for (is = 0; is < QDP_Ns; is++) {
        for (js = 0; js < QDP_Ns; js++) {
            for (ic = 0; ic < arg->nc; ic++) {
                for (jc = 0; jc < arg->nc; jc++, dst++) {
                    QLA_Complex *z = &QLA_elem_P(*src, ic, is, jc, is);
                    QLA_c_eq_r_plus_ir(*dst, QLA_real(*z), QLA_imag(*z));
                }
            }
        }
    }
}

static void
Qs(double_P_get)(char *buf, size_t index, int count, void *v_arg)
{
    Qs(QDPCCArgs) *arg = v_arg;
#if QNc == 'N'
    typedef QLA_DN_DiracPropagator(arg->nc, Ptype);
#else
    typedef Qx(QLA_D,_DiracPropagator) Ptype;
#endif
    Ptype *P = arg->P;
    Ptype *src = &P[index];
    QLA_D_Complex *dst = (void *)buf;
    int ic, jc, is, js;

    if (count != 1)
        luaL_error(arg->L, "qcd.qdpcc.write_prop() count != 1");

    for (is = 0; is < QDP_Ns; is++) {
        for (js = 0; js < QDP_Ns; js++) {
            for (ic = 0; ic < arg->nc; ic++) {
                for (jc = 0; jc < arg->nc; jc++, dst++) {
                    QLA_Complex *z = &QLA_elem_P(*src, ic, is, jc, is);
                    QLA_c_eq_r_plus_ir(*dst, QLA_real(*z), QLA_imag(*z));
                }
            }
        }
    }
}

static int
Qs(cc_write_P_)(lua_State *L, int nc)
{
#if QNc == 'N'
    typedef QLA_DN_DiracPropagator(nc, Ptype);
    typedef QLA_DN_DiracPropagator(nc, Dtype);
    typedef QLA_FN_DiracPropagator(nc, Ftype);
#else
    typedef Qx(QLA_D,_DiracPropagator) Ptype;
    typedef Qx(QLA_D,_DiracPropagator) Dtype;
    typedef Qx(QLA_F,_DiracPropagator) Ftype;
#endif
    typedef void Getter(char *buf, size_t index, int count, void *arg);
/*********************
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
********************/
    const char        prec        = qlua_qio_file_precision(L, 1);
    const char       *fname       = luaL_checkstring(L, 2);
    const char       *file_xml    = luaL_checkstring(L, 3);
    Qs(mLatDirProp)  *prop        = Qs(qlua_checkLatDirProp)(L, 4, NULL, nc);
    mLattice         *S           = qlua_ObjLattice(L, 4);
    int               Sidx        = lua_gettop(L);
    const char       *record_xml  = luaL_checkstring(L, 5);
    int               volfmt      = qlua_qio_volume_format(L, 6, Sidx);
    QIO_Writer       *writer      = NULL;
    Getter           *getter      = NULL;
    char             *file_prec   = NULL;
    int               typesize    = 0;
    int               datacount   = 0;
    int               datum_size  = 0;
    int               word_size   = 0;
    char             *datatype    = NULL;
    QIO_String       *f_xml       = NULL;
    QIO_String       *r_xml       = NULL;
    QIO_RecordInfo   *record_info = NULL;
    Qs(QDPCCArgs) data;

    CALL_QDP(L);
    switch (prec) {
    case 'F':
        file_prec = "F";
        getter = Qs(float_P_get);
        datatype = "QDP_F" Qcolors "_DiracPropagator";
        word_size = sizeof (QLA_F_Real);
        datum_size = sizeof (Ftype);
        break;
    case 'D':
        file_prec = "D";
        getter = Qs(double_P_get);
        datatype = "QDP_D" Qcolors "_DiracPropagator";
        word_size = sizeof (QLA_D_Real);
        datum_size = sizeof (Dtype);
        break;
    default:
        return luaL_error(L, "unsupported output precision");
    }
    datacount = 1;
    typesize = datum_size;

    f_xml = QIO_string_create();
    QIO_string_set(f_xml, file_xml);
    r_xml = QIO_string_create();
    QIO_string_set(r_xml, record_xml);
    writer = qlua_qio_std_writer(S, fname, f_xml, volfmt);
    if (writer == 0)
        return luaL_error(L, "qcd.qdpcc.write_prop(): open error");
    record_info = QIO_create_record_info(QIO_FIELD, NULL, NULL, S->rank,
                                         datatype, file_prec, nc, QDP_Ns,
                                         typesize, datacount);
    data.nc = nc;
    data.P = Qx(QDP_D,_expose_P)(prop->ptr);
    data.L = L;
    if (QIO_write(writer, record_info, r_xml, getter,
                  datum_size, word_size, &data))
        return luaL_error(L, "qcd.qdpcc.write_prop(): record writing error");
        
    if (QIO_close_write(writer))
        return luaL_error(L, "qcd.qdpcc.write_prop(): closing error");

    Qx(QDP_D,_reset_P)(prop->ptr);
    QIO_destroy_record_info(record_info);
    QIO_string_destroy(f_xml);
    QIO_string_destroy(r_xml);

    return 0;
}

#if 0 /* XXX */



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
#endif /* XXX */

#undef Qs
#undef Qx
#undef QC
#undef QNc
#undef Qcolors
