static const char Qs(SeqColMatName)[] = "qcd.ColorMatrix" Qcolors;
static const char Qs(mtnSeqColMat)[] = "qcd.mtColorMatrix" Qcolors;

static int
Qs(q_m_fmt)(lua_State *L)
{
    char fmt[72];
    Qs(mSeqColMat) *b = Qs(qlua_checkSeqColMat)(L, 1, -1);

    sprintf(fmt, "qcd.ColorMatrix%d(%p)", QC(b), b->ptr);
    lua_pushstring(L, fmt);

    return 1;
}

static int
Qs(q_m_gc)(lua_State *L)
{
    Qs(mSeqColMat) *b = Qs(qlua_checkSeqColMat)(L, 1, -1);

    qlua_free(L, b->ptr);
    b->ptr = 0;

    return 0;
}

static int
Qs(q_m_get)(lua_State *L)
{
    switch (qlua_qtype(L, 2)) {
    case qTable: {
        Qs(mSeqColMat) *V = Qs(qlua_checkSeqColMat)(L, 1, -1);
        int b = qlua_checkrightindex(L, 2, QC(V));
        int a = qlua_leftindex(L, 2, QC(V));

        if (a == -1) {
            Qs(mSeqColVec) *r = Qs(qlua_newSeqColVec)(L, QC(V));
                
            Qx(QLA_D,_V_eq_colorvec_M)(QNC(QC(V)) r->ptr, V->ptr, b);
        } else {
            QLA_D_Complex *r = qlua_newComplex(L);
                
            Qx(QLA_D,_C_eq_elem_M)(QNC(QC(V)) r, V->ptr, a, b);
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
Qs(q_m_put)(lua_State *L)
{
    Qs(mSeqColMat) *V = Qs(qlua_checkSeqColMat)(L, 1, -1);
    int b = qlua_checkrightindex(L, 2, QC(V));
    int a = qlua_leftindex(L, 2, QC(V));

    if (a == -1) {
        Qs(mSeqColVec) *z = Qs(qlua_checkSeqColVec)(L, 3, QC(V));
        Qx(QLA_D,_M_eq_colorvec_V)(QNC(QC(V)) V->ptr, z->ptr, b);
    } else {
        QLA_D_Complex *z = qlua_checkComplex(L, 3);

        Qx(QLA_D,_M_eq_elem_C)(QNC(QC(V)) V->ptr, z, a, b);
    }

    return 0;
}

static int
Qs(q_m_norm2_)(lua_State *L)
{
    Qs(mSeqColMat) *a = Qs(qlua_checkSeqColMat)(L, 1, -1);
    QLA_D_Real n;

    Qx(QLA_D,_r_eq_norm2_M)(QNC(QC(a)) &n, a->ptr);
    lua_pushnumber(L, n);
    
    return 1;
}

static int
Qs(q_m_conj)(lua_State *L)
{
    Qs(mSeqColMat) *a = Qs(qlua_checkSeqColMat)(L, 1, -1);
    Qs(mSeqColMat) *r = Qs(qlua_newSeqColMat)(L, QC(a));

    Qx(QLA_D,_M_eq_conj_M)(QNC(QC(a)) r->ptr, a->ptr);

    return 1;
}

static int
Qs(q_m_trans)(lua_State *L)
{
    Qs(mSeqColMat) *a = Qs(qlua_checkSeqColMat)(L, 1, -1);
    Qs(mSeqColMat) *r = Qs(qlua_newSeqColMat)(L, QC(a));

    Qx(QLA_D,_M_eq_transpose_M)(QNC(QC(a)) r->ptr, a->ptr);

    return 1;
}

static int
Qs(q_m_adjoin)(lua_State *L)
{
    Qs(mSeqColMat) *a = Qs(qlua_checkSeqColMat)(L, 1, -1);
    Qs(mSeqColMat) *r = Qs(qlua_newSeqColMat)(L, QC(a));

    Qx(QLA_D,_M_eq_Ma)(QNC(QC(a)) r->ptr, a->ptr);

    return 1;
}

static int
Qs(q_m_trace)(lua_State *L)
{
    Qs(mSeqColMat) *a = Qs(qlua_checkSeqColMat)(L, 1, -1);
    QLA_D_Complex *r = qlua_newComplex(L);

    Qx(QLA_D,_C_eq_trace_M)(QNC(QC(a)) r, a->ptr);

    return 1;
}

static void
#if QNc == 'N'
Qs(XLA_M_eq_exp_M)(int nc,
                   QLA_DN_ColorMatrix(nc, (*r)),
                   QLA_DN_ColorMatrix(nc, (*a)))
#else
Qs(XLA_M_eq_exp_M)(Qx(QLA_D,_ColorMatrix) *r,
                   Qx(QLA_D,_ColorMatrix) *a)
#endif
{
#if QNc == 'N'
    typedef QLA_DN_ColorMatrix(nc, Mtype);
#else
    int nc = QC(a);
    typedef Qx(QLA_D,_ColorMatrix) Mtype;
#endif
    static const int n = 15; /* terms in the Taylor series */
    QLA_D_Complex cone;
    Mtype mone;
    Mtype axs;
    Mtype v;
    int i, j, k;
    double max_el;
    double s;
    unsigned int pw;

    for (max_el = 0, i = 0; i < nc; i++) {
        for (j = 0; j < nc; j++) {
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
    pw = (unsigned int)(ceil(max_el * 2 * nc));
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

static int
Qs(q_m_exp)(lua_State *L)
{
    Qs(mSeqColMat) *a = Qs(qlua_checkSeqColMat)(L, 1, -1);
    Qs(mSeqColMat) *r = Qs(qlua_newSeqColMat)(L, QC(a));

    Qs(XLA_M_eq_exp_M)(QNC(QC(a)) r->ptr, a->ptr);

    return 1;
}

#if (QNc == '2') || (QNc == '3')

/* M_eq_proj_M is adopted from Sergey's qdp-work */
static void
Qs(XLA_sun_proj_step)(Qx(QLA_D,_ColorMatrix) *u, Qx(QLA_D,_ColorMatrix) *w)
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
Qs(XLA_reunit)(Qx(QLA_D,_ColorMatrix) *u)
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
Qs(XLA_reunit)(Qx(QLA_D,_ColorMatrix) *u)
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
/* define XLA_reunit(u) for other number of colors here */
#endif

static void
Qs(XLA_M_eq_proj_M)(Qx(QLA_D,_ColorMatrix) *r,
                    Qx(QLA_D,_ColorMatrix) *a,
                    double BlkAccu,
                    int cnt)
{
    double conver = 1.0;
    QLA_Real new_tr, old_tr;
    Qx(QLA_D,_ColorMatrix) *w = a;

    Qx(QLA_D,_M_eq_M)(r, w);
    Qx(QLA_D,_r_eq_re_M_dot_M)(&new_tr, r, w);
    new_tr = new_tr / QC(xxx);

    while (conver > BlkAccu  && cnt-- > 0) {
        old_tr = new_tr;
        Qs(XLA_sun_proj_step)(r, w);
        Qs(XLA_reunit)(r);
        Qx(QLA_D,_r_eq_re_M_dot_M)(&new_tr, r, w);
        new_tr = new_tr / QC(xxx);
        conver = fabs((new_tr - old_tr) / old_tr);
    }
}

static int
Qs(q_m_proj)(lua_State *L)
{
    Qs(mSeqColMat) *a = Qs(qlua_checkSeqColMat)(L, 1, -1);
    double BlkAccu = luaL_checknumber(L, 2);
    int BlkMax = luaL_checkint(L, 3);
    Qs(mSeqColMat) *r = Qs(qlua_newSeqColMat)(L, QC(a));

    Qs(XLA_M_eq_proj_M)(r->ptr, a->ptr, BlkAccu, BlkMax);

    return 1;
}
#endif

static int
Qs(q_m_set)(lua_State *L)
{
    Qs(mSeqColMat) *r = Qs(qlua_checkSeqColMat)(L, 1, -1);
    Qs(mSeqColMat) *a = Qs(qlua_checkSeqColMat)(L, 2, QC(r));

    Qx(QLA_D,_M_eq_M)(QNC(QC(r)) r->ptr, a->ptr);
    lua_pop(L, 2);

    return 0;
}

static int
Qs(q_m_dot_)(lua_State *L)
{
    Qs(mSeqColMat) *a = Qs(qlua_checkSeqColMat)(L, 1, -1);
    Qs(mSeqColMat) *b = Qs(qlua_checkSeqColMat)(L, 2, QC(a));
    QLA_D_Complex *s = qlua_newComplex(L);

    Qx(QLA_D,_C_eq_M_dot_M)(QNC(QC(a)) s, a->ptr, b->ptr);

    return 1;
}

static int
Qs(q_m_add_m_)(lua_State *L)
{
    Qs(mSeqColMat) *a = Qs(qlua_checkSeqColMat)(L, 1, -1);
    Qs(mSeqColMat) *b = Qs(qlua_checkSeqColMat)(L, 2, QC(a));
    Qs(mSeqColMat) *c = Qs(qlua_newSeqColMat)(L, QC(a));

    Qx(QLA_D,_M_eq_M_plus_M)(QNC(QC(a)) c->ptr, a->ptr, b->ptr);

    return 1;
}

static int
Qs(q_m_sub_m_)(lua_State *L)
{
    Qs(mSeqColMat) *a = Qs(qlua_checkSeqColMat)(L, 1, -1);
    Qs(mSeqColMat) *b = Qs(qlua_checkSeqColMat)(L, 2, QC(a));
    Qs(mSeqColMat) *c = Qs(qlua_newSeqColMat)(L, QC(a));

    Qx(QLA_D,_M_eq_M_minus_M)(QNC(QC(a)) c->ptr, a->ptr, b->ptr);

    return 1;
}

static int
Qs(q_m_neg)(lua_State *L)
{
    Qs(mSeqColMat) *a = Qs(qlua_checkSeqColMat)(L, 1, -1);
    Qs(mSeqColMat) *r = Qs(qlua_newSeqColMat)(L, QC(a));

    Qx(QLA_D,_M_meq_M)(QNC(QC(a))r->ptr, a->ptr);

    return 1;
}

static int
Qs(q_m_mul_m_)(lua_State *L)
{
    Qs(mSeqColMat) *a = Qs(qlua_checkSeqColMat)(L, 1, -1);
    Qs(mSeqColMat) *b = Qs(qlua_checkSeqColMat)(L, 2, QC(a));
    Qs(mSeqColMat) *c = Qs(qlua_newSeqColMat)(L, QC(a));

    Qx(QLA_D,_M_eq_M_times_M)(QNC(QC(a)) c->ptr, a->ptr, b->ptr);
    
    return 1;
}

static int
Qs(q_m_mul_v_)(lua_State *L)
{
    Qs(mSeqColMat) *a = Qs(qlua_checkSeqColMat)(L, 1, -1);
    Qs(mSeqColVec) *b = Qs(qlua_checkSeqColVec)(L, 2, QC(a));
    Qs(mSeqColVec) *c = Qs(qlua_newSeqColVec)(L, QC(a));

    Qx(QLA_D,_V_eq_M_times_V)(QNC(QC(a)) c->ptr, a->ptr, b->ptr);

    return 1;
}

static int
Qs(q_m_mul_r_)(lua_State *L)
{
    Qs(mSeqColMat) *a = Qs(qlua_checkSeqColMat)(L, 1, -1);
    QLA_Real b = luaL_checknumber(L, 2);
    Qs(mSeqColMat) *c = Qs(qlua_newSeqColMat)(L, QC(a));

    Qx(QLA_D,_M_eq_r_times_M)(QNC(QC(a)) c->ptr, &b, a->ptr);

    return 1;
}

static int
Qs(q_r_mul_m_)(lua_State *L)
{
    QLA_Real a = luaL_checknumber(L, 1);
    Qs(mSeqColMat) *b = Qs(qlua_checkSeqColMat)(L, 2, -1);
    Qs(mSeqColMat) *c = Qs(qlua_newSeqColMat)(L, QC(b));

    Qx(QLA_D,_M_eq_r_times_M)(QNC(QC(b)) c->ptr, &a, b->ptr);

    return 1;
}

static int
Qs(q_m_mul_c_)(lua_State *L)
{
    Qs(mSeqColMat) *a = Qs(qlua_checkSeqColMat)(L, 1, -1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    Qs(mSeqColMat) *c = Qs(qlua_newSeqColMat)(L, QC(a));

    Qx(QLA_D,_M_eq_c_times_M)(QNC(QC(a)) c->ptr, b, a->ptr);

    return 1;
}

static int
Qs(q_c_mul_m_)(lua_State *L)
{
    QLA_Complex *a = qlua_checkComplex(L, 1);
    Qs(mSeqColMat) *b = Qs(qlua_checkSeqColMat)(L, 2, -1);
    Qs(mSeqColMat) *c = Qs(qlua_newSeqColMat)(L, QC(b));

    Qx(QLA_D,_M_eq_c_times_M)(QNC(QC(b)) c->ptr, a, b->ptr);

    return 1;
}

static int
Qs(q_m_div_r_)(lua_State *L)
{
    Qs(mSeqColMat) *a = Qs(qlua_checkSeqColMat)(L, 1, -1);
    QLA_D_Real b = 1 / luaL_checknumber(L, 2);
    Qs(mSeqColMat) *c = Qs(qlua_newSeqColMat)(L, QC(a));

    Qx(QLA_D,_M_eq_r_times_M)(QNC(QC(a)) c->ptr, &b, a->ptr);

    return 1;
}

static int
Qs(q_m_div_c_)(lua_State *L)
{
    Qs(mSeqColMat) *a = Qs(qlua_checkSeqColMat)(L, 1, -1);
    QLA_D_Complex *b = qlua_checkComplex(L, 2);
    Qs(mSeqColMat) *c = Qs(qlua_newSeqColMat)(L, QC(a));
    double n = 1 / (QLA_real(*b) * QLA_real(*b) + QLA_imag(*b) * QLA_imag(*b));
    QLA_D_Complex s;

    QLA_real(s) = n * QLA_real(*b);
    QLA_imag(s) = -n * QLA_imag(*b);
    Qx(QLA_D,_M_eq_c_times_M)(QNC(QC(a)) c->ptr, &s, a->ptr);

    return 1;
}

static int
Qs(q_m_colors)(lua_State *L)
{
#if QNc == 'N'
    Qs(mSeqColMat) *a = Qs(qlua_checkSeqColMat)(L, 1, -1);
#else
    Qs(qlua_checkSeqColMat)(L, 1, -1);
#endif
    lua_pushnumber(L, QC(a));

    return 1;
}

static int
Qs(q_m_copy)(lua_State *L)
{
    Qs(mSeqColMat) *a = Qs(qlua_checkSeqColMat)(L, 1, -1);
    Qs(mSeqColMat) *r = Qs(qlua_newSeqColMat)(L, QC(a));

    Qx(QLA_D,_M_eq_M)(QNC(QC(a)) r->ptr, a->ptr);

    return 1;
}

static struct luaL_Reg Qs(mtSeqColMat)[] = {
    { "__tostring",        Qs(q_m_fmt)       },
    { "__gc",              Qs(q_m_gc)        },
    { "__index",           Qs(q_m_get)       },
    { "__newindex",        Qs(q_m_put)       },
    { "__unm",             Qs(q_m_neg)       },
    { "__add",             qlua_add          },
    { "__sub",             qlua_sub          },
    { "__mul",             qlua_mul          },
    { "__div",             qlua_div          },
    { "norm2",             Qs(q_m_norm2_)    },
    { "conj",              Qs(q_m_conj)      },
    { "transpose",         Qs(q_m_trans)     },
    { "adjoin",            Qs(q_m_adjoin)    },
    { "trace",             Qs(q_m_trace)     },
    { "set",               Qs(q_m_set)       },
    { "exp",               Qs(q_m_exp)       },
#if (QNc == '2') || (QNc == '3')
    { "proj",              Qs(q_m_proj)      },
#endif
    { "colors",            Qs(q_m_colors)    },
    { "copy",              Qs(q_m_copy)      },
    /* "a-type" */
    { NULL,                NULL              }
};

Qs(mSeqColMat) *
Qs(qlua_newSeqColMat)(lua_State *L, int nc)
{
#if QNc == 'N'
    typedef QLA_DN_ColorMatrix(nc, DataType);
#else
    typedef Qx(QLA_D,_ColorMatrix) DataType;
#endif
    Qs(mSeqColMat) *hdr = lua_newuserdata(L, sizeof (Qs(mSeqColMat)));
    hdr->ptr = qlua_malloc(L, sizeof (DataType));
#if QNc == 'N'
    hdr->nc = nc;
#endif
    luaL_getmetatable(L, Qs(mtnSeqColMat));
    lua_setmetatable(L, -2);

    return hdr;
}

Qs(mSeqColMat) *
Qs(qlua_newZeroSeqColMat)(lua_State *L, int nc)
{
	Qs(mSeqColMat) *v = Qs(qlua_newSeqColMat)(L, nc);
	Qx(QLA_D,_M_eq_zero)(QNC(nc) v->ptr);
	return v;
}
Qs(mSeqColMat) *
Qs(qlua_checkSeqColMat)(lua_State *L, int idx, int nc)
{
    void *v = qlua_checkLatticeType(L, idx, Qs(qSeqColMat), Qs(SeqColMatName));
    Qs(mSeqColMat) *z = (Qs(mSeqColMat) *)v;

#if QNc == 'N'
    if (nc != -1) {
        if (z->nc != nc)
            luaL_error(L, "Wrong number of colors");
    }
#endif

    return z;
}

static int
Qs(q_seqcolmat_)(lua_State *L, int nc)
{
    switch (lua_gettop(L)) {
    case 0: {
        Qs(qlua_newZeroSeqColMat)(L, nc);
        return 1;
    }
    case 1: {
        switch (qlua_qtype(L, 1)) {
        case qReal: {
            QLA_D_Real x = luaL_checknumber(L, 1);
            Qs(mSeqColMat) *v = Qs(qlua_newSeqColMat)(L, nc);
            QLA_D_Complex z;

            QLA_real(z) = x;
            QLA_imag(z) = 0;
            Qx(QLA_D,_M_eq_c)(QNC(nc) v->ptr, &z);

            return 1;
        }
        case qComplex: {
            QLA_D_Complex *z = qlua_checkComplex(L, 1);
            Qs(mSeqColMat) *v = Qs(qlua_newSeqColMat)(L, nc);

            Qx(QLA_D,_M_eq_c)(QNC(nc) v->ptr, z);

            return 1;
        }
        default:
            break;
        }
        break;
    }
    case 2: {
        switch (qlua_qtype(L, 1)) {
        case qComplex: {
            QLA_D_Complex *z = qlua_checkComplex(L, 1);
            int a = qlua_checkleftindex(L, 2, nc);
            int b = qlua_checkrightindex(L, 2, nc);
            Qs(mSeqColMat) *v = Qs(qlua_newZeroSeqColMat)(L, nc);

            Qx(QLA_D,_M_eq_elem_C)(QNC(nc) v->ptr, z, a, b);

            return 1;
        }
        case Qs(qSeqColVec): {
            Qs(mSeqColVec) *f = Qs(qlua_checkSeqColVec)(L, 1, nc);
            switch (qlua_qtype(L, 2)) {
            case qTable: {
                int b = qlua_checkrightindex(L, 2, nc);
                Qs(mSeqColMat) *v = Qs(qlua_newZeroSeqColMat)(L, nc);
                
                Qx(QLA_D,_M_eq_colorvec_V)(QNC(nc) v->ptr, f->ptr, b);

                return 1;
            }
            case Qs(qSeqColVec): {
                Qs(mSeqColVec) *g = Qs(qlua_checkSeqColVec)(L, 2, nc);
                Qs(mSeqColMat) *v = Qs(qlua_newSeqColMat)(L, nc);

                Qx(QLA_D,_M_eq_V_times_Va)(QNC(nc) v->ptr, f->ptr, g->ptr);

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
    { qlua_add_table, Qs(qSeqColMat),  Qs(qSeqColMat),  Qs(q_m_add_m_) },
    { qlua_sub_table, Qs(qSeqColMat),  Qs(qSeqColMat),  Qs(q_m_sub_m_) },
    { qlua_mul_table, qReal,           Qs(qSeqColMat),  Qs(q_r_mul_m_) },
    { qlua_mul_table, Qs(qSeqColMat),  qReal,           Qs(q_m_mul_r_) },
    { qlua_mul_table, qComplex,        Qs(qSeqColMat),  Qs(q_c_mul_m_) },
    { qlua_mul_table, Qs(qSeqColMat),  qComplex,        Qs(q_m_mul_c_) },
    { qlua_mul_table, Qs(qSeqColMat),  Qs(qSeqColVec),  Qs(q_m_mul_v_) },
    { qlua_mul_table, Qs(qSeqColMat),  Qs(qSeqColMat),  Qs(q_m_mul_m_) },
    { qlua_div_table, Qs(qSeqColMat),  qReal,           Qs(q_m_div_r_) },
    { qlua_div_table, Qs(qSeqColMat),  qComplex,        Qs(q_m_div_c_) },
    { NULL,           qNoType,         qNoType,         NULL           }
};
#undef QNc
#undef Qcolors
#undef Qs
#undef Qx
#undef QC
#undef QNC
