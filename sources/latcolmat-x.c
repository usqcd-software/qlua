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
    Qs(mLatColMat) *b = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);

    Qx(QDP_D,_destroy_M)(b->ptr);
    b->ptr = 0;

    return 0;
}

static int
Qs(q_M_get)(lua_State *L)
{
    switch (qlua_qtype(L, 2)) {
    case qTable: {
        Qs(mLatColMat) *V = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
        mLattice *S = qlua_ObjLattice(L, 1);
        int b = qlua_checkrightindex(L, 2, QC(V));
        int a = qlua_leftindex(L, 2, QC(V));
        int *idx = qlua_latcoord(L, 2, S);
        int Sidx = lua_gettop(L);

        if (idx == NULL) {
            if (a == -1) {
                Qs(mLatColVec) *r = Qs(qlua_newLatColVec)(L, Sidx, QC(V));
                
                CALL_QDP(L);
                Qx(QDP_D,_V_eq_colorvec_M)(r->ptr, V->ptr, b, *S->qss);
            } else {
                mLatComplex *r = qlua_newLatComplex(L, Sidx);
                
                CALL_QDP(L);
                Qx(QDP_D,_C_eq_elem_M)(r->ptr, V->ptr, a, b, *S->qss);
            }
        } else {
            if (a == -1) {
                qlua_free(L, idx);
                return qlua_badindex(L, "ColorMatrix" Qcolors);
            } else {
                QLA_Complex *W = qlua_newComplex(L);
#if QNc == 'N'
                typedef QLA_DN_ColorMatrix(V->nc, Vtype);
#else
                typedef Qx(QLA_D,_ColorMatrix) Vtype;
#endif
                Vtype *locked;
                double zri[2];
                
                qlua_verifylatcoord(L, idx, S);
                CALL_QDP(L);
                locked = Qx(QDP_D,_expose_M)(V->ptr);
                if (QDP_node_number(idx) == QDP_this_node) {
                    QLA_Complex *zz = &Qx(QLA_D,_elem_M)(locked[QDP_index(idx)], a, b);
                    zri[0] = QLA_real(*zz);
                    zri[1] = QLA_imag(*zz);
                } else {
                    zri[0] = 0;
                    zri[1] = 0;
                }
                Qx(QDP_D,_reset_M)(V->ptr);
                QMP_sum_double_array(zri, 2);
                QLA_c_eq_r_plus_ir(*W, zri[0], zri[1]);
            }
        }
        qlua_free(L, idx);
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
    int b = qlua_checkrightindex(L, 2, QC(V));
    int a = qlua_leftindex(L, 2, QC(V));
    int *idx = qlua_latcoord(L, 2, S);

    if (idx == NULL) {
        CALL_QDP(L);
        if (a == -1) {
            Qs(mLatColVec) *z = Qs(qlua_checkLatColVec)(L, 3, S, QC(V));
            Qx(QDP_D,_M_eq_colorvec_V)(V->ptr, z->ptr, b, *S->qss);
        } else {
            mLatComplex *z = qlua_checkLatComplex(L, 3, S);

            Qx(QDP_D,_M_eq_elem_C)(V->ptr, z->ptr, a, b, *S->qss);
        }
    } else {
        if (a == -1) {
            qlua_free(L, idx);
            return qlua_badindex(L, "ColorMatrix" Qcolors);
        } else {
#if QNc == 'N'
            typedef QLA_DN_ColorMatrix(V->nc, Vtype);
#else
            typedef Qx(QLA_D,_ColorMatrix) Vtype;
#endif
            CALL_QDP(L);
            Vtype *locked = Qx(QDP_D,_expose_M)(V->ptr);
            QLA_Complex *zz = &Qx(QLA_D,_elem_M)(locked[QDP_index(idx)], a, b);
            QLA_Complex *z = qlua_checkComplex(L, 3);
            qlua_verifylatcoord(L, idx, S);
            if (QDP_node_number(idx) == QDP_this_node) {
                QLA_c_eq_c(*zz, *z);
            }
            Qx(QDP_D,_reset_M)(V->ptr);
        }
    }
    qlua_free(L, idx);
    return 0;
}

static int
Qs(q_M_norm2_)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_Real n;

    CALL_QDP(L);
    if (S->lss.mask) {
        Qs(mLatColMat) *b = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(a));
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

static struct {
    int nc;
    void *a;
} Qs(Mop_arg); /* YYY global state */

static void
#if QNc == 'N'
Qs(do_Mexp)(int nc, QLA_DN_ColorMatrix(nc, (*r)), int idx)
#else
Qs(do_Mexp)(Qx(QLA_D,_ColorMatrix) *r, int idx)
#endif
{
#if QNc == 'N'
    typedef QLA_DN_ColorMatrix(nc, Mtype);
#else
    typedef Qx(QLA_D,_ColorMatrix) Mtype;
#endif
    static const int n = 15; /* terms in the Taylor series */
    QLA_D_Complex cone;
    Mtype mone;
    Mtype axs;
    Mtype v;
    Mtype *ax = Qs(Mop_arg).a;
    Mtype *a = &ax[idx];
    int i, j, k;
    double max_el;
    double s;
    unsigned int pw;

    for (max_el = 0, i = 0; i < Qs(Mop_arg).nc; i++) {
        for (j = 0; j < Qs(Mop_arg).nc; j++) {
            QLA_D_Complex z;
            double v;
#if QNc == 'N'
            Qx(QLA_D,_C_eq_elem_M)(nc, &z, a, i, j);
#else
            Qx(QLA_D,_C_eq_elem_M)(&z, a, i, j);
#endif
            v = hypot(QLA_real(z), QLA_imag(z));
            if (max_el < v)
                max_el = v;
        }
    }
    pw = (unsigned int)(ceil(max_el * 2 * Qs(Mop_arg).nc));
    QLA_c_eq_r_plus_ir(cone, 1.0, 0.0);
#if QNc == 'N'
    Qx(QLA_D,_M_eq_c)(nc, &mone, &cone);
    s = 1.0 / (n * pw);
    Qx(QLA_D,_M_eq_r_times_M_plus_M)(nc, &v, &s, a, &mone);
    for (k = n; --k;) {
        Qx(QLA_D,_M_eq_M_times_M)(nc, &axs, a, &v);
        s = 1.0 / (k * pw);
        Qx(QLA_D,_M_eq_r_times_M_plus_M)(nc, &v, &s, &axs, &mone);
    }
#else
    Qx(QLA_D,_M_eq_c)(&mone, &cone);
    s = 1.0 / (n * pw);
    Qx(QLA_D,_M_eq_r_times_M_plus_M)(&v, &s, a, &mone);
    for (k = n; --k;) {
        Qx(QLA_D,_M_eq_M_times_M)(&axs, a, &v);
        s = 1.0 / (k * pw);
        Qx(QLA_D,_M_eq_r_times_M_plus_M)(&v, &s, &axs, &mone);
    }
#endif

#if QNc == 'N'
    Qx(QLA_D,_M_eq_M)(nc, r, &mone);
    while (pw) {
        if (pw & 1) {
            Qx(QLA_D,_M_eq_M_times_M)(nc, &axs, r, &v);
            Qx(QLA_D,_M_eq_M)(nc, r, &axs);
        }
        pw >>= 1;
        if (pw > 0) {
            Qx(QLA_D,_M_eq_M_times_M)(nc, &axs, &v, &v);
            Qx(QLA_D,_M_eq_M)(nc, &v, &axs);
        }
    }
#else
    Qx(QLA_D,_M_eq_M)(r, &mone);
    while (pw) {
        if (pw & 1) {
            Qx(QLA_D,_M_eq_M_times_M)(&axs, r, &v);
            Qx(QLA_D,_M_eq_M)(r, &axs);
        }
        pw >>= 1;
        if (pw > 0) {
            Qx(QLA_D,_M_eq_M_times_M)(&axs, &v, &v);
            Qx(QLA_D,_M_eq_M)(&v, &axs);
        }
    }
#endif
}

static void
Qs(X_M_eq_exp_M)(int nc,
                 Qx(QDP_D,_ColorMatrix) *r,
                 Qx(QDP_D,_ColorMatrix) *a,
                 QDP_Subset s)
{
    Qs(Mop_arg).nc = nc;
    Qs(Mop_arg).a = Qx(QDP_D,_expose_M)(a);
    Qx(QDP_D,_M_eq_funci)(r, Qs(do_Mexp), s);
    Qx(QDP_D,_reset_M)(a);
    Qs(Mop_arg).a = 0;
    Qs(Mop_arg).nc = -1;
}

static int
Qs(q_M_exp)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColMat) *r = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qs(X_M_eq_exp_M)(QC(a), r->ptr, a->ptr, *S->qss);

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
    int i1, i2;

    for (i1 = 1; i1 < QC(xxx); i1++) {
        for (i2 = 0; i2 < i1; i2++) {
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
        QLA_Complex *ui0 = &Mu(i, 0);
        t = hypot(t, QLA_real(*ui0));
        t = hypot(t, QLA_imag(*ui0));
    }
    t = 1.0 / t;
    /* rescale the first column */
    for (i = 0; i < QC(xxx); i++) {
        QLA_Complex *ui0 = &Mu(i, 0);

        QLA_real(*ui0) *= t;
        QLA_imag(*ui0) *= t;
    }
    /* compute u[0]* . u[1] */
    QLA_c_eq_r_plus_ir(t0k, 0.0, 0.0);
    for (i = 0; i < QC(xxx); i++)
        QLA_c_peq_ca_times_c(t0k, Mu(i, 0), Mu(i, 1));

    /* make u[1] orthogonal to u[0]: u[1] = u[1] - t0k * u[0] */
    for (i = 0; i < QC(xxx); i++)
        QLA_c_meq_c_times_c(Mu(i, 1), t0k, Mu(i, 0));

    /* normalize u[1] */
    for (t = 0, i = 0; i < QC(xxx); i++) {
        QLA_Complex *ui1 = &Mu(i, 1);
        t = hypot(t, QLA_real(*ui1));
        t = hypot(t, QLA_imag(*ui1));
    }
    t = 1.0 / t;
    for (i = 0; i < QC(xxx); i++) {
        QLA_Complex *ui1 = &Mu(i, 1);

        QLA_real(*ui1) *= t;
        QLA_imag(*ui1) *= t;
    }

    /* u[2] is a vector product of u[0]* and u[1]* */
    QLA_c_eq_ca_times_ca (Mu(0,2), Mu(1,0), Mu(2,1));
    QLA_c_meq_ca_times_ca(Mu(0,2), Mu(2,0), Mu(1,1));
    
    QLA_c_eq_ca_times_ca (Mu(1,2), Mu(2,0), Mu(0,1));
    QLA_c_meq_ca_times_ca(Mu(1,2), Mu(0,0), Mu(2,1));

    QLA_c_eq_ca_times_ca (Mu(2,2), Mu(0,0), Mu(1,1));
    QLA_c_meq_ca_times_ca(Mu(2,2), Mu(1,0), Mu(0,1));
#undef Mu
}
#else
/* define X_reunit(u) for other number of colors here */
#endif

static struct {
    Qx(QLA_D,_ColorMatrix) *a;
    double BlkAccu;
    int BlkMax;
} Qs(Mproj_args);

static void
Qs(do_Mproj)(Qx(QLA_D,_ColorMatrix) *r, int idx)
{
    int cnt = Qs(Mproj_args).BlkMax;
    double BlkAccu = Qs(Mproj_args).BlkAccu;
    double conver = 1.0;
    QLA_Real new_tr, old_tr;
    Qx(QLA_D,_ColorMatrix) *w = &Qs(Mproj_args).a[idx];

    Qx(QLA_D,_M_eq_M)(r, w);
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

static void
Qs(X_M_eq_proj_M)(Qx(QDP_D,_ColorMatrix) *r,
                  Qx(QDP_D,_ColorMatrix) *a,
                  double BlkAccu, int BlkMax, QDP_Subset s)
{
    Qs(Mproj_args).a = Qx(QDP_D,_expose_M)(a);
    Qs(Mproj_args).BlkAccu = BlkAccu;
    Qs(Mproj_args).BlkMax = BlkMax;
    Qx(QDP_D,_M_eq_funci)(r, Qs(do_Mproj), s);
    Qx(QDP_D,_reset_M)(a);
    Qs(Mop_arg).a = 0;
    Qs(Mproj_args).BlkAccu = 0;
    Qs(Mproj_args).BlkMax = 0;
}

static int
Qs(q_M_proj)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    double BlkAccu = luaL_checknumber(L, 2);
    int BlkMax = luaL_checkint(L, 3);
    Qs(mLatColMat) *r = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qs(X_M_eq_proj_M)(r->ptr, a->ptr, BlkAccu, BlkMax, *S->qss);

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
    if (S->lss.mask)
        Qx(QDP_D,_M_eq_M_mask_I)(r->ptr, a->ptr, S->lss.mask, *S->qss);
    else
        Qx(QDP_D,_M_eq_M)(r->ptr, a->ptr, *S->qss);
    lua_pop(L, 2);

    return 1;
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
    Qx(QDP_D,_M_meq_M)(r->ptr, a->ptr, *S->qss);

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

static struct {
    int nc;
    QLA_D_Real *a;
    void *b; /* Qx(QLA_D,_ColorMatrix) */
} Qs(RMmul_args); /* YYY global state */

#if QNc == 'N'
static void
Qs(do_RMmul)(int nc, QLA_DN_ColorMatrix(nc, (*r)), int idx)
{
    QLA_DN_ColorMatrix(nc, (*b)) = Qs(RMmul_args).b;
    Qx(QLA_D,_M_eq_r_times_M)(nc, r, &Qs(RMmul_args).a[idx], &b[idx]);
}
#else
static void
Qs(do_RMmul)(Qx(QLA_D,_ColorMatrix) *r, int idx)
{
    Qx(QLA_D,_ColorMatrix) *b = Qs(RMmul_args).b;
    Qx(QLA_D,_M_eq_r_times_M)(r, &Qs(RMmul_args).a[idx], &b[idx]);
}
#endif

static void
Qs(X_M_eq_R_times_M)(int nc,
                     Qx(QDP_D,_ColorMatrix) *r,
                     QDP_D_Real *a,
                     Qx(QDP_D, _ColorMatrix) *b,
                     QDP_Subset s)
{
    Qs(RMmul_args).nc = nc;
    Qs(RMmul_args).a = QDP_D_expose_R(a);
    Qs(RMmul_args).b = Qx(QDP_D,_expose_M)(b);
    Qx(QDP_D,_M_eq_funci)(r, Qs(do_RMmul), s);
    Qx(QDP_D,_reset_M)(b);
    QDP_D_reset_R(a);
    Qs(RMmul_args).a = 0;
    Qs(RMmul_args).b = 0;
    Qs(RMmul_args).nc = -1;
}

static int
Qs(q_R_mul_M_)(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColMat) *b = Qs(qlua_checkLatColMat)(L, 2, S, -1);
    Qs(mLatColMat) *r = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(b));

    CALL_QDP(L);
    Qs(X_M_eq_R_times_M)(QC(b), r->ptr, a->ptr, b->ptr, *S->qss);

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
    Qs(X_M_eq_R_times_M)(QC(a), r->ptr, b->ptr, a->ptr, *S->qss);

    return 1;
}

static struct {
    int nc;
    QLA_D_Complex *a;
    void *b; /* Qx(QLA_D,_ColorMatrix) */
} Qs(CMmul_args); /* YYY global state */

#if QNc == 'N'
static void
Qs(do_CMmul)(int nc, QLA_DN_ColorMatrix(nc, (*r)), int idx)
{
    QLA_DN_ColorMatrix(nc, (*b)) = Qs(CMmul_args).b;
    Qx(QLA_D,_M_eq_c_times_M)(nc, r, &Qs(CMmul_args).a[idx], &b[idx]);
}
#else
static void
Qs(do_CMmul)(Qx(QLA_D,_ColorMatrix) *r, int idx)
{
    Qx(QLA_D,_ColorMatrix) *b = Qs(CMmul_args).b;
    Qx(QLA_D,_M_eq_c_times_M)(r, &Qs(CMmul_args).a[idx], &b[idx]);
}
#endif

static void
Qs(X_M_eq_C_times_M)(int nc,
                     Qx(QDP_D,_ColorMatrix) *r,
                     QDP_D_Complex *a,
                     Qx(QDP_D,_ColorMatrix) *b,
                     QDP_Subset s)
{
    Qs(CMmul_args).nc = nc;
    Qs(CMmul_args).a = QDP_D_expose_C(a);
    Qs(CMmul_args).b = Qx(QDP_D,_expose_M)(b);
    Qx(QDP_D,_M_eq_funci)(r, Qs(do_CMmul), s);
    Qx(QDP_D,_reset_M)(b);
    QDP_D_reset_C(a);
    Qs(CMmul_args).a = 0;
    Qs(CMmul_args).b = 0;
    Qs(CMmul_args).nc = -1;
}

static int
Qs(q_C_mul_M_)(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColMat) *b = Qs(qlua_checkLatColMat)(L, 2, S, -1);
    Qs(mLatColMat) *r = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(b));

    CALL_QDP(L);
    Qs(X_M_eq_C_times_M)(QC(b), r->ptr, a->ptr, b->ptr, *S->qss);

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
    Qs(X_M_eq_C_times_M)(QC(a), r->ptr, b->ptr, a->ptr, *S->qss);

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
    { "norm2",             Qs(q_M_norm2_)    },
    { "shift",             Qs(q_M_shift)     },
    { "conj",              Qs(q_M_conj)      },
    { "transpose",         Qs(q_M_trans)     },
    { "adjoin",            Qs(q_M_adjoin)    },
    { "trace",             Qs(q_M_trace)     },
    { "set",               Qs(q_M_set)       },
    { "exp",               Qs(q_M_exp)       },
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
#if QNc == 'N'
    Qx(QDP_D,_ColorMatrix) *v = Qx(QDP_D,_create_M)(nc);
#else
    Qx(QDP_D,_ColorMatrix) *v = Qx(QDP_D,_create_M)();
#endif
    Qs(mLatColMat) *hdr;

    if (v == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
#if QNc == 'N'
        v = Qx(QDP_D,_create_M)(nc);
#else
        v = Qx(QDP_D,_create_M)();
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

    return hdr;
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
Qs(q_latcolmat_)(lua_State *L, mLattice *S, int nc, int off)
{
    switch (lua_gettop(L) - off) {
    case 1: {
        Qs(mLatColMat) *v = Qs(qlua_newLatColMat)(L, 1, nc);

        CALL_QDP(L);
        Qx(QDP_D,_M_eq_zero)(v->ptr, *S->qss);
        
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
            Qs(mLatColMat) *v = Qs(qlua_newLatColMat)(L, 1, nc);

            CALL_QDP(L);
            Qx(QDP_D,_M_eq_zero)(v->ptr, *S->qss);
            Qx(QDP_D,_M_eq_elem_C)(v->ptr, z->ptr, a, b, *S->qss);

            return 1;
        }
        case Qs(qLatColVec): {
            Qs(mLatColVec) *f = Qs(qlua_checkLatColVec)(L, 2 + off, S, nc);
            switch (qlua_qtype(L, 3 + off)) {
            case qTable: {
                int b = qlua_checkrightindex(L, 3 + off, nc);
                Qs(mLatColMat) *v = Qs(qlua_newLatColMat)(L, 1, nc);
                
                CALL_QDP(L);
                Qx(QDP_D,_M_eq_zero)(v->ptr, *S->qss);
                Qx(QDP_D,_M_eq_colorvec_V)(v->ptr, f->ptr, b, *S->qss);

                return 1;
            }
            case Qs(qLatColVec): {
                Qs(mLatColVec) *g = Qs(qlua_checkLatColVec)(L, 3 + off, S, nc);
                Qs(mLatColMat) *v = Qs(qlua_newLatColMat)(L, 1, nc);

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

#undef QNc
#undef Qcolors
#undef Qs
#undef Qx
#undef QC
