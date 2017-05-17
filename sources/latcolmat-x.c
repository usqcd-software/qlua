static const char Qs(LatColMatName)[] = "lattice.ColorMatrix" Qcolors;

static int
Qs(q_M_fmt)(lua_State *L)
{
    char fmt[72];
    Qs(mLatColMat) *b = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);

    sprintf(fmt, "QDP:ColorMatrix%d(%p)", QC(b), b->ptr);
    lua_pushstring(L, fmt);

    return 1;
}

static int
Qs(q_M_gc)(lua_State *L)
{
    char qdp_name[72];
    Qs(mLatColMat) *b = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);

    sprintf(qdp_name, "ColorMatrix%d", QC(b));
    qlua_qdp_memuse(L, qdp_name, -1);
    Qx(QDP_D,_destroy_M)(b->ptr);
    b->ptr = 0;

    return 0;
}

static void
#if QNc == 'N'
Qs(unpack_mat)(int nc, double *std, QLA_DN_ColorMatrix(nc, (*src)))
#else
Qs(unpack_mat)(int nc, double *std, Qx(QLA_D,_ColorMatrix) *src)
#endif
{
    int a, b;

    for (a = 0; a < nc; a++) {
        for (b = 0; b < nc; b++) {
            QLA_D_Complex *z = &Qx(QLA_D,_elem_M)(*src, a, b);
            std[2 * (a + nc * b)] = QLA_real(*z);
            std[2 * (a + nc * b) + 1] = QLA_imag(*z);
        }
    }
}

static void
#if QNc == 'N'
Qs(pack_mat)(int nc, QLA_DN_ColorMatrix(nc, (*dst)), const double *std)
#else
Qs(pack_mat)(int nc, Qx(QLA_D,_ColorMatrix) *dst, const double *std)
#endif
{
    QLA_D_Complex z;
    int a, b;

    for (a = 0; a < nc; a++) {
        for (b = 0; b < nc; b++) {
            int idx = 2 * (a + nc * b);
            QLA_c_eq_r_plus_ir(z, std[idx], std[idx + 1]);
            QLA_c_eq_c(Qx(QLA_D,_elem_M)(*dst, a, b), z);
        }
    }
}

static int
Qs(q_M_get)(lua_State *L)
{
    switch (qlua_qtype(L, 2)) {
    case qTable: {
        Qs(mLatColMat) *V = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
        mLattice *S = qlua_ObjLattice(L, 1);
        int *idx = qlua_latcoord(L, 2, S);
        int Sidx = lua_gettop(L);
#if QNc == 'N'
        typedef QLA_DN_ColorMatrix(QC(V), Vtype);
#else
        typedef Qx(QLA_D,_ColorMatrix) Vtype;
#endif

        if (idx == NULL) {
            int b = qlua_checkrightindex(L, 2, QC(V));
            int a = qlua_leftindex(L, 2, QC(V));
            if (a == -1) {
                /* M[{b=expr}] => V */
                Qs(mLatColVec) *r = Qs(qlua_newLatColVec)(L, Sidx, QC(V));
                
                CALL_QDP(L);
                Qx(QDP_D,_V_eq_colorvec_M)(r->ptr, V->ptr, b, *S->qss);
            } else {
                /* M[{a=expr,b=expr}] => C */
                mLatComplex *r = qlua_newLatComplex(L, Sidx);
                
                CALL_QDP(L);
                Qx(QDP_D,_C_eq_elem_M)(r->ptr, V->ptr, a, b, *S->qss);
            }
        } else {
            qlua_verifylatcoord(L, idx, S);
            int b = qlua_rightindex(L, 2, QC(V));
            int a = qlua_leftindex(L, 2, QC(V));
            CALL_QDP(L);
            Vtype *locked = Qx(QDP_D,_expose_M)(V->ptr);
            int site_node = QDP_node_number_L(S->lat, idx);

            if ((a == -1) && (b == -1)) {
                Qs(mSeqColMat) *m = Qs(qlua_newSeqColMat)(L, QC(V));
                double *m_std = qlua_malloc(L, 2 * QC(V) * QC(V) * sizeof (double));

                if (site_node == QDP_this_node) {
                    int la = QDP_index_L(S->lat, idx);
                    Qs(unpack_mat)(QC(V), m_std, &locked[la]);
                }
                XMP_dist_double_array(site_node, 2 * QC(V) * QC(V), m_std);
                Qs(pack_mat)(QC(V), m->ptr, m_std);
                                qlua_free(L, m_std);
            } else if ((a == -1) || (b == -1)) {
                qlua_free(L, idx);
                return qlua_badindex(L, "ColorMatrix" Qcolors);
            } else {
                /* M[{x,y,z...,a=expr,b=expr}] => c */
                QLA_Complex *W = qlua_newComplex(L);
                double zri[2];
                
                if (site_node == QDP_this_node) {
                    int la = QDP_index_L(S->lat, idx);
                    QLA_Complex *z = &Qx(QLA_D,_elem_M)(locked[la], a, b);
                    zri[0] = QLA_real(*z);
                    zri[1] = QLA_imag(*z);
                }
                XMP_dist_double_array(site_node, 2, zri);
                QLA_c_eq_r_plus_ir(*W, zri[0], zri[1]);
            }
            Qx(QDP_D,_reset_M)(V->ptr);
            qlua_free(L, idx);
        }
        return 1;
    }
    case qString:
        return qlua_selflookup(L, 1, luaL_checkstring(L, 2));
    default:
        break;
    }
    return qlua_badindex(L, "ColorMatrix" Qcolors);
}

static int
Qs(q_M_put)(lua_State *L)
{
    Qs(mLatColMat) *V = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    int *idx = qlua_latcoord(L, 2, S);

    if (idx == NULL) {
        int b = qlua_checkrightindex(L, 2, QC(V));
        int a = qlua_leftindex(L, 2, QC(V));

        CALL_QDP(L);
        if (a == -1) {
            /* M[{b=expr}] = V */
            Qs(mLatColVec) *z = Qs(qlua_checkLatColVec)(L, 3, S, QC(V));
            Qx(QDP_D,_M_eq_colorvec_V)(V->ptr, z->ptr, b, *S->qss);
        } else {
            /* M[{a=expr,b=expr}] = C */
            mLatComplex *z = qlua_checkLatComplex(L, 3, S);

            Qx(QDP_D,_M_eq_elem_C)(V->ptr, z->ptr, a, b, *S->qss);
        }
    } else {
        qlua_verifylatcoord(L, idx, S);
        int b = qlua_rightindex(L, 2, QC(V));
        int a = qlua_leftindex(L, 2, QC(V));
#if QNc == 'N'
        typedef QLA_DN_ColorMatrix(QC(V), Vtype);
#else
        typedef Qx(QLA_D,_ColorMatrix) Vtype;
#endif
        CALL_QDP(L);
        Vtype *locked = Qx(QDP_D,_expose_M)(V->ptr);
        int site_node = QDP_node_number_L(S->lat, idx);

        if ((a == -1) && (b == -1)) {
            /* M[{x,y,z...}] = m */
            Qs(mSeqColMat) *v = Qs(qlua_checkSeqColMat)(L, 3, QC(V));
            if (site_node == QDP_this_node) {
                int la = QDP_index_L(S->lat, idx);
                Qx(QLA_D,_M_eq_M)(QNC(QC(V)) &locked[la], v->ptr);
            }
        } else if ((a == -1) || (b == -1)) {
            qlua_free(L, idx);
            return qlua_badindex(L, "ColorMatrix" Qcolors);
        } else {
            /* M[{x,y,z...,a=expr,b=expr}] = c */
            QLA_D_Complex *z = qlua_checkComplex(L, 3);
            if (site_node == QDP_this_node) {
                int la = QDP_index_L(S->lat, idx);
                QLA_D_Complex *w = &Qx(QLA_D,_elem_M)(locked[la], a, b);
                QLA_c_eq_c(*w, *z);
            }
        }
        Qx(QDP_D,_reset_M)(V->ptr);
        qlua_free(L, idx);
    }
    return 0;
}

static int
Qs(q_M_sum)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    int argc = lua_gettop(L);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    int nc = QC(a);
    
    switch (argc) {
    case 1: {
        Qs(mSeqColMat) *s = Qs(qlua_newSeqColMat)(L, nc);

        CALL_QDP(L);
        if (S->lss.mask) {
            Qs(mLatColMat) *b = Qs(qlua_newZeroLatColMat)(L, Sidx, nc);
            Qx(QDP_D,_M_eq_M_mask_I)(b->ptr, a->ptr, S->lss.mask, *S->qss);
            Qx(QDP_D,_m_eq_sum_M)(s->ptr, b->ptr, *S->qss);
            lua_pop(L, 1);
        } else {
            Qx(QDP_D,_m_eq_sum_M)(s->ptr, a->ptr, *S->qss);
        }
        return 1;
    }
    case 2: {
#if QNc == 'N'
        typedef QLA_DN_ColorMatrix(nc, Vtype);
#else
        typedef Qx(QLA_D,_ColorMatrix) Vtype;
#endif
        mLatMulti *m = qlua_checkLatMulti(L, 2, S);
        int size = m->size;
        QLA_Int *ii = m->idx;
        int sites = QDP_sites_on_node_L(S->lat);
        Vtype **vv = qlua_malloc(L, size * sizeof(Vtype *));
        int i, k, ca, cb;

        lua_createtable(L, size, 0);
        for (i = 0; i < size; i++) {
            Qs(mSeqColMat) *vi = Qs(qlua_newZeroSeqColMat)(L, QC(a));
            vv[i] = vi->ptr;
            lua_rawseti(L, -2, i + 1); /* [sic] lua index */
        }
        CALL_QDP(L);
        Vtype *xx = Qx(QDP_D,_expose_M)(a->ptr);
        
        for (k = 0; k < sites; k++, xx++, ii++) {
            int t = *ii;
            if ((t < 0) || (t >= size))
                continue;
            Qx(QLA_D,_M_peq_M)(QNC(nc) vv[t], xx);
        }
        Qx(QDP_D,_reset_M)(a->ptr);
        QLA_D_Real *rr = qlua_malloc(L, 2 * size * nc * nc * sizeof (QLA_D_Real));
        for (i = 0; i < size; i++) {
            for (ca = 0; ca < nc; ca++) {
                for (cb = 0; cb < nc; cb++) {
                    QLA_D_Complex *z = &Qx(QLA_D,_elem_M)(*vv[i], ca, cb);
                    int ab = 2 * (ca + nc * (cb + nc * i));
                    rr[ab] = QLA_real(*z);
                    rr[ab + 1] = QLA_imag(*z);
                }
            }
        }
        QMP_sum_double_array(rr, 2 * size * nc * nc);
        for (i = 0; i < size; i++) {
            for (ca = 0; ca < nc; ca++) {
                for (cb = 0; cb < nc; cb++) {
                    QLA_D_Complex z;
                    int ab =  2 * (ca + nc * (cb + nc * i));
                    QLA_c_eq_r_plus_ir(z, rr[ab], rr[ab + 1]);
                    QLA_c_eq_c(Qx(QLA_D,_elem_M)(*vv[i], ca, cb), z);
                }
            }
        }
                qlua_free(L, rr);
                qlua_free(L, vv);
        return 1;
    }
    }
    return luaL_error(L, "bad arguments for ColorMatrix:sum()");
}

static int
Qs(q_M_norm2_)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_D_Real n = 0;

    CALL_QDP(L);
    if (S->lss.mask) {
        Qs(mLatColMat) *b = Qs(qlua_newZeroLatColMat)(L, lua_gettop(L), QC(a));
        Qx(QDP_D,_M_eq_M_mask_I)(b->ptr, a->ptr, S->lss.mask, *S->qss);
        Qx(QDP_D,_r_eq_norm2_M)(&n, b->ptr, *S->qss);
    } else {
        Qx(QDP_D,_r_eq_norm2_M)(&n, a->ptr, *S->qss);
    }

    lua_pushnumber(L, n);
    
    return 1;
}

static int
Qs(q_M_shift)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    QDP_Shift shift = qlua_checkShift(L, 2, S);
    QDP_ShiftDir dir = qlua_checkShiftDir(L, 3);
    Qs(mLatColMat) *r = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_M_eq_sM)(r->ptr, a->ptr, shift, dir, *S->qss);

    return 1;
}

static int
Qs(q_M_conj)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColMat) *r = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_M_eq_conj_M)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_M_trans)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColMat) *r = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_M_eq_transpose_M)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_M_adjoin)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColMat) *r = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_M_eq_Ma)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_M_trace)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    Qx(QDP_D,_C_eq_trace_M)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_M_exp)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColMat) *r = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_M_eq_exp_M)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_M_det)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    Qx(QDP_D,_C_eq_det_M)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_M_log)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColMat) *r = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_M_eq_log_M)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_M_sqrt)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColMat) *r = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_M_eq_sqrt_M)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_M_inverse)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColMat) *r = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_M_eq_inverse_M)(r->ptr, a->ptr, *S->qss);

    return 1;
}

#if (QNc == '2') || (QNc == '3')

/* M_eq_proj_M is adopted from Sergey's qdp-work */
static void
Qs(X_sun_proj_step)(Qx(QLA_D,_ColorMatrix) *u, Qx(QLA_D,_ColorMatrix) *w)
{
    Qx(QLA_D,_ColorMatrix) v;
    Qx(QLA_D,_ColorMatrix) tmp;
    double r0, r1, r2, r3;
    int di, i1, i2;

    for (di = 1 ; di < QC(xxx) ; di++) {
        for (i1 = 0 ; i1 < QC(xxx) - di ; i1++) {
            i2 = i1 + di;
            int j;
            double r_l;
            double a0, a1, a2, a3;

            Qx(QLA_D,_M_eq_M_times_Ma)(&v, u, w);
            /* extract su(2) part at (i,j) */
            r0 = (QLA_real(Qx(QLA_D,_elem_M)(v,i1,i1)) +
                  QLA_real(Qx(QLA_D,_elem_M)(v,i2,i2)));
            r1 = (QLA_imag(Qx(QLA_D,_elem_M)(v,i1,i2)) +
                  QLA_imag(Qx(QLA_D,_elem_M)(v,i2,i1)));
            r2 = (QLA_real(Qx(QLA_D,_elem_M)(v,i1,i2)) -
                  QLA_real(Qx(QLA_D,_elem_M)(v,i2,i1)));
            r3 = (QLA_imag(Qx(QLA_D,_elem_M)(v,i1,i1)) -
                  QLA_imag(Qx(QLA_D,_elem_M)(v,i2,i2)));
            r_l = hypot(hypot(r0, r1), hypot(r2, r3));
            if (r_l > COLOR_EPS) {
                r_l = 1.0 / r_l;
                a0 = r0 * r_l;
                a1 = -r1 * r_l;
                a2 = -r2 * r_l;
                a3 = -r3 * r_l;
            } else {
                a0 = 1.0; a1 = a2 = a3 = 0.0;
            }
            /* put a's back into v */
            Qx(QLA_D,_M_eq_zero)(&v);
            for (j = 0; j < QC(xxx); j++)
                QLA_c_eq_r_plus_ir(Qx(QLA_D,_elem_M)(v,j,j), 1.0, 0.0);
            QLA_c_eq_r_plus_ir(Qx(QLA_D,_elem_M)(v,i1,i1),  a0,  a3);
            QLA_c_eq_r_plus_ir(Qx(QLA_D,_elem_M)(v,i1,i2),  a2,  a1);
            QLA_c_eq_r_plus_ir(Qx(QLA_D,_elem_M)(v,i2,i1), -a2,  a1);
            QLA_c_eq_r_plus_ir(Qx(QLA_D,_elem_M)(v,i2,i2),  a0, -a3);
            
            Qx(QLA_D,_M_eq_M_times_M)(&tmp, &v, u);
            Qx(QLA_D,_M_eq_M)(u, &tmp);
        }
    }
}

#if QNc == '2'
static void
Qs(X_reunit)(Qx(QLA_D,_ColorMatrix) *u)
{
    QLA_Complex *u00 = &Qx(QLA_D,_elem_M)(*u, 0, 0);
    QLA_Complex *u10 = &Qx(QLA_D,_elem_M)(*u, 1, 0);
    double t = 1.0 / hypot(hypot(QLA_real(*u00), QLA_imag(*u00)),
                           hypot(QLA_real(*u10), QLA_imag(*u10)));
    double u00r = QLA_real(*u00) * t;
    double u00i = QLA_imag(*u00) * t;
    double u10r = QLA_real(*u10) * t;
    double u10i = QLA_imag(*u10) * t;

    QLA_c_eq_r_plus_ir(Qx(QLA_D,_elem_M)(*u, 0, 0),  u00r,  u00i);
    QLA_c_eq_r_plus_ir(Qx(QLA_D,_elem_M)(*u, 1, 0),  u10r,  u10i);
    QLA_c_eq_r_plus_ir(Qx(QLA_D,_elem_M)(*u, 0, 1), -u10r,  u10i);
    QLA_c_eq_r_plus_ir(Qx(QLA_D,_elem_M)(*u, 1, 1),  u00r, -u00i);
}
#elif QNc == '3'
static void
Qs(X_reunit)(Qx(QLA_D,_ColorMatrix) *u)
{
    QLA_Complex t0k;
    double t;
    int i;

#define Mu(i,j)   Qx(QLA_D,_elem_M)(*u, i, j)
    /* norm of the first column */
    for (t = 0, i = 0; i < QC(xxx); i++) {
        QLA_Complex *ui0 = &Mu(0, i);
        t = hypot(t, QLA_real(*ui0));
        t = hypot(t, QLA_imag(*ui0));
    }
    t = 1.0 / t;
    /* rescale the first column */
    for (i = 0; i < QC(xxx); i++) {
        QLA_Complex *ui0 = &Mu(0, i);

        QLA_real(*ui0) *= t;
        QLA_imag(*ui0) *= t;
    }
    /* compute u[0]* . u[1] */
    QLA_c_eq_r_plus_ir(t0k, 0.0, 0.0);
    for (i = 0; i < QC(xxx); i++)
        QLA_c_peq_ca_times_c(t0k, Mu(0, i), Mu(1, i));

    /* make u[1] orthogonal to u[0]: u[1] = u[1] - t0k * u[0] */
    for (i = 0; i < QC(xxx); i++)
        QLA_c_meq_c_times_c(Mu(1, i), t0k, Mu(0, i));

    /* normalize u[1] */
    for (t = 0, i = 0; i < QC(xxx); i++) {
        QLA_Complex *ui1 = &Mu(1, i);
        t = hypot(t, QLA_real(*ui1));
        t = hypot(t, QLA_imag(*ui1));
    }
    t = 1.0 / t;
    for (i = 0; i < QC(xxx); i++) {
        QLA_Complex *ui1 = &Mu(1, i);

        QLA_real(*ui1) *= t;
        QLA_imag(*ui1) *= t;
    }

    /* u[2] is a vector product of u[0]* and u[1]* */
    QLA_c_eq_ca_times_ca (Mu(2,0), Mu(0,1), Mu(1,2));
    QLA_c_meq_ca_times_ca(Mu(2,0), Mu(0,2), Mu(1,1));
                                                  
    QLA_c_eq_ca_times_ca (Mu(2,1), Mu(0,2), Mu(1,0));
    QLA_c_meq_ca_times_ca(Mu(2,1), Mu(0,0), Mu(1,2));
                                                  
    QLA_c_eq_ca_times_ca (Mu(2,2), Mu(0,0), Mu(1,1));
    QLA_c_meq_ca_times_ca(Mu(2,2), Mu(0,1), Mu(1,0));
#undef Mu
}
#else
/* define X_reunit(u) for other number of colors here */
#endif

typedef struct {
    Qx(QLA_D,_ColorMatrix) *a;
    double BlkAccu;
    int BlkMax;
} Qs(Mproj_arg);

static void
Qs(do_Mproj)(Qx(QLA_D,_ColorMatrix) *r, int idx, void *env)
{
    Qs(Mproj_arg) *arg = env;
    int cnt = arg->BlkMax;
    double BlkAccu = arg->BlkAccu;
    double conver = 1.0;
    QLA_Real new_tr, old_tr;
    Qx(QLA_D,_ColorMatrix) *w = &arg->a[idx];

    Qs(X_reunit)(r); /* r is initalized to initial value for iter.proj */
    Qx(QLA_D,_r_eq_re_M_dot_M)(&new_tr, r, w);
    new_tr = new_tr / QC(xxx);

    while (conver > BlkAccu  && cnt-- > 0) {
        old_tr = new_tr;
        Qs(X_sun_proj_step)(r, w);
        Qs(X_reunit)(r);
        Qx(QLA_D,_r_eq_re_M_dot_M)(&new_tr, r, w);
        new_tr = new_tr / QC(xxx);
        conver = fabs((new_tr - old_tr) / old_tr);
    }
}

static int
Qs(q_M_proj)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    double BlkAccu = luaL_checknumber(L, 2);
    int BlkMax = luaL_checkint(L, 3);
    Qs(mLatColMat) *r0= (3 < lua_gettop(L) 
                         ? Qs(qlua_checkLatColMat)(L, 4, NULL, -1)
                         : NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColMat) *r = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(a));
    if (NULL != r0) /* specified initial val. for iterative projection*/
        Qx(QDP_D,_M_eq_M)(r->ptr, r0->ptr, *S->qss);
    else            /* start with the matrix itself as initial; 
                       this is likely incorrect as the matrix is not SU(N) */
        Qx(QDP_D,_M_eq_M)(r->ptr, a->ptr, *S->qss);
    
    Qs(Mproj_arg) arg;

    arg.a = Qx(QDP_D,_expose_M)(a->ptr);
    arg.BlkAccu = BlkAccu;
    arg.BlkMax = BlkMax;
    CALL_QDP(L);
    Qx(QDP_D,_M_eq_funcia)(r->ptr, Qs(do_Mproj), &arg, *S->qss);
    Qx(QDP_D,_reset_M)(a->ptr);

    return 1;
}
#endif

static int
Qs(q_M_set)(lua_State *L)
{
    Qs(mLatColMat) *r = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 2, S, QC(r));

    CALL_QDP(L);
    if (S->lss.mask) {
        Qx(QDP_D,_M_eq_M_mask_I)(r->ptr, a->ptr, S->lss.mask, *S->qss);
    } else {
        Qx(QDP_D,_M_eq_M)(r->ptr, a->ptr, *S->qss);
    }
    lua_pop(L, 2);

    return 0;
}

static int
Qs(q_M_dot_)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColMat) *b = Qs(qlua_checkLatColMat)(L, 2, S, QC(a));
    mLatComplex *s = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    Qx(QDP_D,_C_eq_M_dot_M)(s->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
Qs(q_M_add_M_)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColMat) *b = Qs(qlua_checkLatColMat)(L, 2, S, QC(a));
    Qs(mLatColMat) *c = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_M_eq_M_plus_M)(c->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
Qs(q_M_sub_M_)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColMat) *b = Qs(qlua_checkLatColMat)(L, 2, S, QC(a));
    Qs(mLatColMat) *c = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_M_eq_M_minus_M)(c->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
Qs(q_M_neg)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColMat) *r = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_M_eqm_M)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_M_mul_M_)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColMat) *b = Qs(qlua_checkLatColMat)(L, 2, S, QC(a));
    Qs(mLatColMat) *c = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_M_eq_M_times_M)(c->ptr, a->ptr, b->ptr, *S->qss);
    
    return 1;
}

static int
Qs(q_M_mul_V_)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColVec) *b = Qs(qlua_checkLatColVec)(L, 2, S, QC(a));
    Qs(mLatColVec) *c = Qs(qlua_newLatColVec)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_V_eq_M_times_V)(c->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
Qs(q_M_mul_r_)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_Real b = luaL_checknumber(L, 2);
    Qs(mLatColMat) *c = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_M_eq_r_times_M)(c->ptr, &b, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_r_mul_M_)(lua_State *L)
{
    QLA_Real a = luaL_checknumber(L, 1);
    Qs(mLatColMat) *b = Qs(qlua_checkLatColMat)(L, 2, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 2);
    Qs(mLatColMat) *c = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(b));

    CALL_QDP(L);
    Qx(QDP_D,_M_eq_r_times_M)(c->ptr, &a, b->ptr, *S->qss);

    return 1;
}

static int
Qs(q_M_mul_c_)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    Qs(mLatColMat) *c = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_M_eq_c_times_M)(c->ptr, b, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_c_mul_M_)(lua_State *L)
{
    QLA_Complex *a = qlua_checkComplex(L, 1);
    Qs(mLatColMat) *b = Qs(qlua_checkLatColMat)(L, 2, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 2);
    Qs(mLatColMat) *c = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(b));

    CALL_QDP(L);
    Qx(QDP_D,_M_eq_c_times_M)(c->ptr, a, b->ptr, *S->qss);

    return 1;
}

static int
Qs(q_R_mul_M_)(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColMat) *b = Qs(qlua_checkLatColMat)(L, 2, S, -1);
    Qs(mLatColMat) *r = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(b));

    CALL_QDP(L);
    Qx(QDP_D,_M_eq_R_times_M)(r->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
Qs(q_M_mul_R_)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, S);
    Qs(mLatColMat) *r = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_M_eq_R_times_M)(r->ptr, b->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_C_mul_M_)(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColMat) *b = Qs(qlua_checkLatColMat)(L, 2, S, -1);
    Qs(mLatColMat) *r = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(b));

    CALL_QDP(L);
    Qx(QDP_D,_M_eq_C_times_M)(r->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
Qs(q_M_mul_C_)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2, S);
    Qs(mLatColMat) *r = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_M_eq_C_times_M)(r->ptr, b->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_M_div_r_)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_D_Real b = 1 / luaL_checknumber(L, 2);
    Qs(mLatColMat) *c = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_M_eq_r_times_M)(c->ptr, &b, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_M_div_c_)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_D_Complex *b = qlua_checkComplex(L, 2);
    Qs(mLatColMat) *c = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(a));
    double n = 1 / (QLA_real(*b) * QLA_real(*b) + QLA_imag(*b) * QLA_imag(*b));
    QLA_D_Complex s;

    CALL_QDP(L);
    QLA_real(s) = n * QLA_real(*b);
    QLA_imag(s) = -n * QLA_imag(*b);
    Qx(QDP_D,_M_eq_c_times_M)(c->ptr, &s, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_M_colors)(lua_State *L)
{
#if QNc == 'N'
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
#else
    Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
#endif
    lua_pushnumber(L, QC(a));

    return 1;
}

static int
Qs(q_M_copy)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColMat) *r = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_M_eq_M)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static struct luaL_Reg Qs(mtLatColMat)[] = {
    { "__tostring",        Qs(q_M_fmt)       },
    { "__gc",              Qs(q_M_gc)        },
    { "__index",           Qs(q_M_get)       },
    { "__newindex",        Qs(q_M_put)       },
    { "__unm",             Qs(q_M_neg)       },
    { "__add",             qlua_add          },
    { "__sub",             qlua_sub          },
    { "__mul",             qlua_mul          },
    { "__div",             qlua_div          },
    { "sum",               Qs(q_M_sum)       },
    { "norm2",             Qs(q_M_norm2_)    },
    { "shift",             Qs(q_M_shift)     },
    { "conj",              Qs(q_M_conj)      },
    { "transpose",         Qs(q_M_trans)     },
    { "adjoin",            Qs(q_M_adjoin)    },
    { "trace",             Qs(q_M_trace)     },
    { "set",               Qs(q_M_set)       },
    { "det",               Qs(q_M_det)       },
    { "exp",               Qs(q_M_exp)       },
    { "log",               Qs(q_M_log)       },
    { "sqrt",              Qs(q_M_sqrt)      },
    { "inverse",           Qs(q_M_inverse)   },
#if (QNc == '2') || (QNc == '3')
    { "proj",              Qs(q_M_proj)      },
#endif
    { "colors",            Qs(q_M_colors)    },
    { "copy",              Qs(q_M_copy)      },
    /* "lattice" */
    /* "a-type" */
    { NULL,                NULL              }
};

Qs(mLatColMat) *
Qs(qlua_newLatColMat)(lua_State *L, int Sidx, int nc)
{
    mLattice *S = qlua_checkLattice(L, Sidx);
    char qdp_name[72];
#if QNc == 'N'
    Qx(QDP_D,_ColorMatrix) *v = Qx(QDP_D,_create_M_L)(nc, S->lat);
#else
    Qx(QDP_D,_ColorMatrix) *v = Qx(QDP_D,_create_M_L)(S->lat);
#endif
    Qs(mLatColMat) *hdr;

    if (v == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
#if QNc == 'N'
        v = Qx(QDP_D,_create_M_L)(nc, S->lat);
#else
        v = Qx(QDP_D,_create_M_L)(S->lat);
#endif
        if (v == 0)
            luaL_error(L, "not enough memory (QDP_ColorMatrix" Qcolors ")");
    }
    hdr = lua_newuserdata(L, sizeof (Qs(mLatColMat)));
    hdr->ptr = v;
#if QNc == 'N'
    hdr->nc = nc;
#endif
    qlua_createLatticeTable(L, Sidx, Qs(mtLatColMat), Qs(qLatColMat),
                            Qs(LatColMatName));
    lua_setmetatable(L, -2);
    sprintf(qdp_name, "ColorMatrix%d", QC(hdr));
    qlua_qdp_memuse(L, qdp_name, 1);

    return hdr;
}

Qs(mLatColMat) *
Qs(qlua_newZeroLatColMat)(lua_State *L, int Sidx, int nc)
{
        Qs(mLatColMat) *v = Qs(qlua_newLatColMat)(L, Sidx, nc);
        mLattice *S = qlua_checkLattice(L, Sidx);
        Qx(QDP_D,_M_eq_zero)(v->ptr, S->all);
        return v;
}

Qs(mLatColMat) *
Qs(qlua_checkLatColMat)(lua_State *L, int idx, mLattice *S, int nc)
{
    void *v = qlua_checkLatticeType(L, idx, Qs(qLatColMat), Qs(LatColMatName));
    Qs(mLatColMat) *z;

    if (S) {
        mLattice *S1 = qlua_ObjLattice(L, idx);
        if (S1->id != S->id)
            luaL_error(L, "%s on a wrong lattice", Qs(LatColMatName));
        lua_pop(L, 1);
    }
    z = (Qs(mLatColMat) *)v;
#if QNc == 'N'
    if (nc != -1) {
        if (z->nc != nc)
            luaL_error(L, "Wrong number of colors");
    }
#endif

    return z;
}

static int
Qs(q_latcolmat_seq_)(lua_State *L, mLattice *S, int nc)
{
    Qs(mSeqColMat) *m = Qs(qlua_checkSeqColMat)(L, 2, nc);
    Qs(mLatColMat) *M = Qs(qlua_newLatColMat)(L, 1, nc);

    CALL_QDP(L);
    Qx(QDP_D, _M_eq_m)(M->ptr, m->ptr, *S->qss);

    return 1;
}

static int
Qs(q_latcolmat_lat_)(lua_State *L, mLattice *S, int nc)
{
    Qs(mLatColMat) *m = Qs(qlua_checkLatColMat)(L, 2, S, nc);
    Qs(mLatColMat) *M = Qs(qlua_newLatColMat)(L, 1, nc);

    CALL_QDP(L);
    Qx(QDP_D, _M_eq_M)(M->ptr, m->ptr, *S->qss);

    return 1;
}

static int
Qs(q_latcolmat_)(lua_State *L, mLattice *S, int nc, int off)
{
    switch (lua_gettop(L) - off) {
    case 1: {
        Qs(qlua_newZeroLatColMat)(L, 1, nc);
        return 1;
    }
    case 2: {
        switch (qlua_qtype(L, 2 + off)) {
        case qReal: {
            QLA_D_Real x = luaL_checknumber(L, 2 + off);
            Qs(mLatColMat) *v = Qs(qlua_newLatColMat)(L, 1, nc);
            QLA_D_Complex z;

            CALL_QDP(L);
            QLA_real(z) = x;
            QLA_imag(z) = 0;
            Qx(QDP_D,_M_eq_c)(v->ptr, &z, *S->qss);

            return 1;
        }
        case qComplex: {
            QLA_D_Complex *z = qlua_checkComplex(L, 2 + off);
            Qs(mLatColMat) *v = Qs(qlua_newLatColMat)(L, 1, nc);

            CALL_QDP(L);
            Qx(QDP_D,_M_eq_c)(v->ptr, z, *S->qss);

            return 1;
        }
        default:
            break;
        }
        break;
    }
    case 3: {
        switch (qlua_qtype(L, 2 + off)) {
        case qLatComplex: {
            mLatComplex *z = qlua_checkLatComplex(L, 2 + off, S);
            int a = qlua_checkleftindex(L, 3 + off, nc);
            int b = qlua_checkrightindex(L, 3 + off, nc);
            Qs(mLatColMat) *v = Qs(qlua_newZeroLatColMat)(L, 1, nc);

            CALL_QDP(L);
            Qx(QDP_D,_M_eq_elem_C)(v->ptr, z->ptr, a, b, *S->qss);

            return 1;
        }
        case Qs(qLatColVec): {
            Qs(mLatColVec) *f = Qs(qlua_checkLatColVec)(L, 2 + off, S, nc);
            switch (qlua_qtype(L, 3 + off)) {
            case qTable: {
                int b = qlua_checkrightindex(L, 3 + off, nc);
                Qs(mLatColMat) *v = Qs(qlua_newZeroLatColMat)(L, 1, nc);
                
                CALL_QDP(L);
                Qx(QDP_D,_M_eq_colorvec_V)(v->ptr, f->ptr, b, *S->qss);

                return 1;
            }
            case Qs(qLatColVec): {
                Qs(mLatColVec) *g = Qs(qlua_checkLatColVec)(L, 3 + off, S, nc);
                Qs(mLatColMat) *v = Qs(qlua_newZeroLatColMat)(L, 1, nc);

                CALL_QDP(L);
                Qx(QDP_D,_M_eq_V_times_Va)(v->ptr, f->ptr, g->ptr, *S->qss);

                return 1;
            }
            default:
                break;
            }
            break;
        }
        default:
            break;
        }
        break;
    }
    }
    return qlua_badconstr(L, "ColorMatrix" Qcolors);
}

static const QLUA_Op2 Qs(ops)[] = {
    { qlua_add_table, Qs(qLatColMat),  Qs(qLatColMat),  Qs(q_M_add_M_) },
    { qlua_sub_table, Qs(qLatColMat),  Qs(qLatColMat),  Qs(q_M_sub_M_) },
    { qlua_mul_table, qReal,           Qs(qLatColMat),  Qs(q_r_mul_M_) },
    { qlua_mul_table, Qs(qLatColMat),  qReal,           Qs(q_M_mul_r_) },
    { qlua_mul_table, qComplex,        Qs(qLatColMat),  Qs(q_c_mul_M_) },
    { qlua_mul_table, Qs(qLatColMat),  qComplex,        Qs(q_M_mul_c_) },
    { qlua_mul_table, qLatReal,        Qs(qLatColMat),  Qs(q_R_mul_M_) },
    { qlua_mul_table, Qs(qLatColMat),  qLatReal,        Qs(q_M_mul_R_) },
    { qlua_mul_table, qLatComplex,     Qs(qLatColMat),  Qs(q_C_mul_M_) },
    { qlua_mul_table, Qs(qLatColMat),  qLatComplex,     Qs(q_M_mul_C_) },
    { qlua_mul_table, Qs(qLatColMat),  Qs(qLatColVec),  Qs(q_M_mul_V_) },
    { qlua_mul_table, Qs(qLatColMat),  Qs(qLatColMat),  Qs(q_M_mul_M_) },
    { qlua_div_table, Qs(qLatColMat),  qReal,           Qs(q_M_div_r_) },
    { qlua_div_table, Qs(qLatColMat),  qComplex,        Qs(q_M_div_c_) },
    { NULL,           qNoType,         qNoType,         NULL           }
};

#undef QNc
#undef Qcolors
#undef Qs
#undef Qx
#undef QC
#undef QNC
