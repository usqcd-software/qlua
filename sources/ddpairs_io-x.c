typedef struct {
    int         nc;
    lua_State  *L;
    void       *P;     /* QLA_DiracPropagator   *P; */
    int        js;
    int        jc;
} Qs(USQCDArgs);

static void
Qs(float_DD_put)(char *buf, size_t index, int count, void *v_arg)
{
    Qs(USQCDArgs) *arg = v_arg;
#if QNc == 'N'
    typedef QLA_DN_DiracPropagator(arg->nc, Ptype);
#else
    typedef Qx(QLA_D,_DiracPropagator) Ptype;
#endif
    Ptype *P = arg->P;
    Ptype *dst = &P[index];
    QLA_F_Complex *src = (void *)buf;
    int js = arg->js;
    int jc = arg->jc;
    int is, ic;
    
    if (count != 1)
        luaL_error(arg->L, "qcd.ddpairs.read(): count != 1");
    
    for (is = 0; is < QDP_Ns; is++) {
        for (ic = 0; ic < arg->nc; ic++, src++) {
            QLA_c_eq_r_plus_ir(QLA_elem_P(*dst, ic, is, jc, js),
                               QLA_real(*src), QLA_imag(*src));
        }
    }
}

static void
Qs(double_DD_put)(char *buf, size_t index, int count, void *v_arg)
{
    Qs(USQCDArgs) *arg = v_arg;
#if QNc == 'N'
    typedef QLA_DN_DiracPropagator(arg->nc, Ptype);
#else
    typedef Qx(QLA_D,_DiracPropagator) Ptype;
#endif
    Ptype *P = arg->P;
    QLA_D_Complex *src = (void *)buf;
    Ptype *dst = &P[index];
    int js = arg->js;
    int jc = arg->jc;
    int is, ic;
    
    if (count != 1)
        luaL_error(arg->L, "qcd.ddpairs.read(): count != 1");
    
    for (is = 0; is < QDP_Ns; is++) {
        for (ic = 0; ic < arg->nc; ic++, src++) {
            QLA_c_eq_r_plus_ir(QLA_elem_P(*dst, ic, is, jc, js),
                               QLA_real(*src), QLA_imag(*src));
        }
    }
}

static int
Qs(dd_read_comp)(int nc,
                 lua_State *L,
                 int js, int jc,
#if QNc == 'N'
                 QLA_DN_DiracPropagator(nc, (*T)),
#else
                 Qx(QLA_D,_DiracPropagator) *T,
#endif
                 QIO_Reader *reader)
{
#if QNc == 'N'
    typedef QLA_DN_DiracPropagator(nc, Ptype);
    typedef QLA_DN_DiracFermion(nc, Dtype);
    typedef QLA_FN_DiracFermion(nc, Ftype);
#else
    typedef Qx(QLA_D,_DiracPropagator) Ptype;
    typedef Qx(QLA_D,_DiracFermion) Dtype;
    typedef Qx(QLA_F,_DiracFermion) Ftype;
#endif
    typedef void Putter(char *buf, size_t index, int count, void *arg);
    QIO_String      *record_xml = NULL;
    QIO_RecordInfo  *record_info = NULL;
    int              word_size = 0;
    int              datum_size = 0;
    Putter          *putter = NULL;
    int              prec;
    Qs(USQCDArgs)    data;

    record_info = QIO_create_record_info(0, NULL, NULL, 0, "", "", 0, 0, 0, 0);
    record_xml = QIO_string_create();
    if (QIO_read_record_info(reader, record_info, record_xml))
        return luaL_error(L, "qcd.ddpairs.read() read info error");
    prec = *QIO_get_precision(record_info);
    switch (prec) {
    case 'F':
        putter = Qs(float_DD_put);
        word_size = sizeof (QLA_F_Real);
        datum_size = sizeof (Ftype);
        break;
    case 'D':
        putter = Qs(double_DD_put);
        word_size = sizeof (QLA_D_Real);
        datum_size = sizeof (Dtype);
        break;
    default:
        return luaL_error(L, "unsupported file precision");
    }

    data.nc = nc;
    data.P = T;
    data.js = js;
    data.jc = jc;
    data.L = L;
    if (QIO_read(reader, record_info, record_xml, putter,
                 datum_size, word_size, &data))
        return luaL_error(L, "qcd.ddpairs.read() read error");

    QIO_string_destroy(record_xml);
    QIO_destroy_record_info(record_info);
    return 0;
}

static int
Qs(dd_read)(lua_State *L, int off, int nc)
{
#if QNc == 'N'
    typedef QLA_DN_DiracPropagator(nc, Ptype);
#else
    typedef Qx(QLA_D,_DiracPropagator) Ptype;
#endif
    mLattice          *S        = qlua_checkLattice(L, off + 1);
    const char        *fname    = luaL_checkstring(L, off + 2);
    Qs(mLatDirProp)   *src      = Qs(qlua_newLatDirProp)(L, off + 1, nc);
    Qs(mLatDirProp)   *prop     = Qs(qlua_newLatDirProp)(L, off + 1, nc);
    QIO_String        *file_xml = NULL;
    QIO_Reader        *reader   = NULL;
    int js, jc;
    Ptype *xP;
    Ptype *xS;
    
    CALL_QDP(L);
    file_xml = QIO_string_create();
    reader = qlua_qio_std_reader(S, fname, file_xml);
    if (reader == 0)
        return luaL_error(L, "qcd.ddpairs.read() open error");
    lua_pushstring(L, QIO_string_ptr(file_xml));
    QIO_string_destroy(file_xml);
    xP = Qx(QDP_D,_expose_P)(prop->ptr);
    xS = Qx(QDP_D,_expose_P)(src->ptr);
    for (js = 0; js < QDP_Ns; js++) {
        for (jc = 0; jc < nc; jc++) {
            Qs(dd_read_comp)(nc, L, js, jc, xS, reader);
            Qs(dd_read_comp)(nc, L, js, jc, xP, reader);
        }
    }
    Qx(QDP_D,_reset_P)(prop->ptr);
    Qx(QDP_D,_reset_P)(src->ptr);
    QIO_close_read(reader);
    return 3;
}

static void
Qs(float_DD_get)(char *buf, size_t index, int count, void *v_arg)
{
    Qs(USQCDArgs) *arg = v_arg;
#if QNc == 'N'
    typedef QLA_DN_DiracPropagator(arg->nc, Ptype);
#else
    typedef Qx(QLA_D,_DiracPropagator) Ptype;
#endif
    Ptype *P = arg->P;
    Ptype *src = &P[index];
    QLA_F_Complex *dst = (void *)buf;
    int jc = arg->jc;
    int js = arg->js;
    int ic, is;

    for (is = 0; is < QDP_Ns; is++) {
        for (ic = 0; ic < arg->nc; ic++, dst++) {
            QLA_D_Complex *z = &QLA_elem_P(*src, ic, is, jc, js);
            QLA_c_eq_r_plus_ir(*dst, QLA_real(*z), QLA_imag(*z));
        }
    }
}

static void
Qs(double_DD_get)(char *buf, size_t index, int count, void *v_arg)
{
    Qs(USQCDArgs) *arg = v_arg;
#if QNc == 'N'
    typedef QLA_DN_DiracPropagator(arg->nc, Ptype);
#else
    typedef Qx(QLA_D,_DiracPropagator) Ptype;
#endif
    Ptype *P = arg->P;
    Ptype *src = &P[index];
    QLA_D_Complex *dst = (void *)buf;
    int jc = arg->jc;
    int js = arg->js;
    int ic, is;

    for (is = 0; is < QDP_Ns; is++) {
        for (ic = 0; ic < arg->nc; ic++, dst++) {
            QLA_D_Complex *z = &QLA_elem_P(*src, ic, is, jc, js);
            QLA_c_eq_r_plus_ir(*dst, QLA_real(*z), QLA_imag(*z));
        }
    }
}

static int
Qs(dd_write)(lua_State *L, int nc)
{
#if QNc == 'N'
    typedef QLA_DN_DiracPropagator(nc, Ptype);
    typedef QLA_DN_DiracFermion(nc, Dtype);
    typedef QLA_FN_DiracFermion(nc, Ftype);
#else
    typedef Qx(QLA_D,_DiracPropagator) Ptype;
    typedef Qx(QLA_D,_DiracFermion) Dtype;
    typedef Qx(QLA_F,_DiracFermion) Ftype;
#endif
    typedef void Getter(char *buf, size_t index, int count, void *arg);
    const char        prec        = qlua_qio_file_precision(L, 1);
    const char       *fname       = luaL_checkstring(L, 2);
    const char       *file_xml    = luaL_checkstring(L, 3);
    Qs(mLatDirProp)  *source      = Qs(qlua_checkLatDirProp)(L, 4, NULL, nc);
    mLattice         *S           = qlua_ObjLattice(L, 4);
    int               Sidx        = lua_gettop(L);
    const char       *source_xml  = luaL_checkstring(L, 5);
    int               time_slice  = luaL_checkint(L, 6);
    Qs(mLatDirProp)  *prop        = Qs(qlua_checkLatDirProp)(L, 7, S, nc);
    const char       *prop_xml    = luaL_checkstring(L, 8);
    int               volfmt      = qlua_qio_volume_format(L, 9, Sidx);
    QIO_Writer       *writer      = NULL;
    char             *file_prec   = NULL;
    Getter           *getter      = NULL;
    int               typesize    = 0;
    int               datacount   = 0;
    int               datum_size  = 0;
    int               word_size   = 0;
    char             *datatype    = NULL;
    QIO_String       *f_xml       = NULL;
    QIO_String       *s_xml       = NULL;
    QIO_String       *p_xml       = NULL;
    int               lo[S->rank];
    int               hi[S->rank];
    Qs(USQCDArgs)     dataP;
    Qs(USQCDArgs)     dataS;
    int               i, jc, js;

    CALL_QDP(L);
    for (i = 0; i < S->rank; i++) {
        lo[i] = 0;
        hi[i] = S->dim[i] - 1;
    }
    if ((time_slice >= 0) && (time_slice < S->dim[S->rank - 1]))
        lo[S->rank - 1] = hi[S->rank - 1] = time_slice;
    switch (prec) {
    case 'F':
        file_prec = "F";
        getter = Qs(float_DD_get);
        datatype = "QDP_F" Qcolors "_DiracFermion";
        word_size = sizeof (QLA_F_Real);
        datum_size = sizeof (Ftype);
        break;
    case 'D':
        file_prec = "D";
        getter = Qs(double_DD_get);
        datatype = "QDP_D" Qcolors "_DiracFermion";
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
    s_xml = QIO_string_create();
    QIO_string_set(s_xml, source_xml);
    p_xml = QIO_string_create();
    QIO_string_set(p_xml, prop_xml);
    writer = qlua_qio_std_writer(S, fname, f_xml, volfmt);
    if (writer == 0)
        return luaL_error(L, "qcd.ddpairs.write(): open error");
    dataP.nc = nc;
    dataP.P = Qx(QDP_D,_expose_P)(prop->ptr);
    dataP.L = L;
    dataS.nc = nc;
    dataS.P = Qx(QDP_D,_expose_P)(source->ptr);
    dataS.L = L;
    for (js = 0; js < QDP_Ns; js++) {
        for (jc = 0; jc < nc; jc++) {
            QIO_RecordInfo *rec_info;
            int rec_cl; /* This is how qio.h defines it */
            int *rec_lo, *rec_hi;

            dataS.js = js;
            dataS.jc = jc;
            if ((time_slice >= 0) && (time_slice < S->dim[S->rank - 1])) {
                rec_cl = QIO_HYPER;
                rec_lo = lo;
                rec_hi = hi;
            } else {
                rec_cl = QIO_FIELD;
                rec_lo = NULL;
                rec_hi = NULL;
            }
            rec_info = QIO_create_record_info(rec_cl, rec_lo, rec_hi, S->rank,
                                              datatype, file_prec, nc, QDP_Ns,
                                              typesize, datacount);
            if (QIO_write(writer, rec_info, s_xml, getter,
                          datum_size, word_size, &dataS))
                return luaL_error(L, "qcd.ddpairs.write(): source write error");
            QIO_destroy_record_info(rec_info);
            
            dataP.js = js;
            dataP.jc = jc;
            rec_info = QIO_create_record_info(QIO_FIELD, NULL, NULL, S->rank,
                                              datatype, file_prec, nc, QDP_Ns,
                                              typesize, datacount);
            if (QIO_write(writer, rec_info, p_xml, getter,
                          datum_size, word_size, &dataP))
                return luaL_error(L, "qcd.ddpairs.write(): prop write error");
            QIO_destroy_record_info(rec_info);
        }
    }

    Qx(QDP_D,_reset_P)(source->ptr);
    Qx(QDP_D,_reset_P)(prop->ptr);
    if (QIO_close_write(writer))
        return luaL_error(L, "qcd.ddpairs.write(): closing error");
    QIO_string_destroy(p_xml);
    QIO_string_destroy(s_xml);
    QIO_string_destroy(f_xml);

    return 0;
}

#undef Qs
#undef Qx
#undef QC
#undef QNc
#undef Qcolors
