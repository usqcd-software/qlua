
const char Qs(LatDirPropName)[] = "lattice.DiracPropagator" Qcolors;

static int
Qs(q_P_fmt)(lua_State *L)
{
    char fmt[72];
    Qs(mLatDirProp) *b = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);

    sprintf(fmt, "QDP:DiracPropagator%d(%p)", QC(b), b->ptr);
    lua_pushstring(L, fmt);

    return 1;
}

static int
Qs(q_P_gc)(lua_State *L)
{
    Qs(mLatDirProp) *b = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);

#if QNc == 'N'
    Qx(QDP_D,_destroy_P)(QC(b), b->ptr);
#else
    Qx(QDP_D,_destroy_P)(b->ptr);
#endif
    b->ptr = 0;

    return 0;
}

static int
Qs(q_P_get)(lua_State *L)
{
    switch (qlua_qtype(L, 2)) {
    case qTable: {
        Qs(mLatDirProp) *V = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
        mLattice *S = qlua_ObjLattice(L, 1);
        int d = qlua_checkdiracindex(L, 2);
        int c = qlua_checkcolorindex(L, 2, QC(V));
        Qs(mLatDirFerm) *r = Qs(qlua_newLatDirFerm)(L, lua_gettop(L), QC(V));
                
        CALL_QDP(L);
#if QNc == 'N'
        Qx(QDP_D,_D_eq_diracvec_P)(QC(V), r->ptr, V->ptr, c, d, *S->qss);
#else
        Qx(QDP_D,_D_eq_diracvec_P)(r->ptr, V->ptr, c, d, *S->qss);
#endif

        return 1;
    }
    case qString:
        return qlua_selflookup(L, 1, luaL_checkstring(L, 2));
    default:
        break;
    }
    return qlua_badindex(L, "DiracPropagator" Qcolors);
}

static int
Qs(q_P_put)(lua_State *L)
{
    Qs(mLatDirProp) *V = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    int d = qlua_checkdiracindex(L, 2);
    int c = qlua_checkcolorindex(L, 2, QC(V));
    Qs(mLatDirFerm) *r = Qs(qlua_checkLatDirFerm)(L, 3, S, QC(V));
            
    CALL_QDP(L);
#if QNc == 'N'
    Qx(QDP_D,_P_eq_diracvec_D)(QC(V), V->ptr, r->ptr, c, d, *S->qss);
#else
    Qx(QDP_D,_P_eq_diracvec_D)(V->ptr, r->ptr, c, d, *S->qss);
#endif

    return 0;
}

static int
Qs(q_P_norm2_)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_D_Real n;

    CALL_QDP(L);
#if QNc == 'N'
    Qx(QDP_D,_r_eq_norm2_P)(QC(a), &n, a->ptr, *S->qss);
#else
    Qx(QDP_D,_r_eq_norm2_P)(&n, a->ptr, *S->qss);
#endif
    lua_pushnumber(L, n);
    
    return 1;
}

static int
Qs(q_P_shift)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    QDP_Shift shift = qlua_checkShift(L, 2, S);
    QDP_ShiftDir dir = qlua_checkShiftDir(L, 3);
    Qs(mLatDirProp) *r = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
#if QNc == 'N'
    Qx(QDP_D,_P_eq_sP)(QC(a), r->ptr, a->ptr, shift, dir, *S->qss);
#else
    Qx(QDP_D,_P_eq_sP)(r->ptr, a->ptr, shift, dir, *S->qss);
#endif

    return 1;
}

static int
Qs(q_P_conj)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirProp) *r = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
#if QNc == 'N'
    Qx(QDP_D,_P_eq_conj_P)(QC(a), r->ptr, a->ptr, *S->qss);
#else
    Qx(QDP_D,_P_eq_conj_P)(r->ptr, a->ptr, *S->qss);
#endif

    return 1;
}

static int
Qs(q_P_trans)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirProp) *r = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
#if QNc == 'N'
    Qx(QDP_D,_P_eq_transpose_P)(QC(a), r->ptr, a->ptr, *S->qss);
#else
    Qx(QDP_D,_P_eq_transpose_P)(r->ptr, a->ptr, *S->qss);
#endif

    return 1;
}

static int
Qs(q_P_adjoin)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirProp) *r = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
#if QNc == 'N'
    Qx(QDP_D,_P_eq_Pa)(QC(a), r->ptr, a->ptr, *S->qss);
#else
    Qx(QDP_D,_P_eq_Pa)(r->ptr, a->ptr, *S->qss);
#endif

    return 1;
}

static int
Qs(q_P_spintrace)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColMat) *r = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
#if QNc == 'N'
    Qx(QDP_D,_M_eq_spintrace_P)(QC(a), r->ptr, a->ptr, *S->qss);
#else
    Qx(QDP_D,_M_eq_spintrace_P)(r->ptr, a->ptr, *S->qss);
#endif

    return 1;
}

static struct {
    int nc;
    Qx(QLA_D,_DiracPropagator) *a;
} Qs(Pst_args); /* YYY global state */

static void
Qs(do_Pst)(Qx(QLA_D,_DiracPropagator) *r, int idx)
{
    Qx(QLA_D,_DiracPropagator) *a = &Qs(Pst_args).a[idx];
    int nc = Qs(Pst_args).nc;
    int is, ic, js, jc;

    for (is = 0; is < QDP_Ns; is++) {
        for (js = 0; js < QDP_Ns; js++) {
            for (ic = 0; ic < nc; ic++)  {
                for (jc = 0; jc < nc; jc++) {
                    QLA_c_eq_c(Qx(QLA_D,_elem_P)(*r,ic,is,jc,js),
                               Qx(QLA_D,_elem_P)(*a,ic,js,jc,is));
                }
            }
        }
    }
}

static int
Qs(q_P_spintranspose)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirProp) *r = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qs(Pst_args).nc = QC(a);
#if QNc == 'N'
    Qs(Pst_args).a = Qx(QDP_D,_expose_P)(QC(a), a->ptr);
    Qx(QDP_D,_P_eq_funci)(QC(a), r->ptr, Qs(do_Pst), *S->qss);
    Qx(QDP_D,_reset_P)(QC(a), a->ptr);
#else
    Qs(Pst_args).a = Qx(QDP_D,_expose_P)(a->ptr);
    Qx(QDP_D,_P_eq_funci)(r->ptr, Qs(do_Pst), *S->qss);
    Qx(QDP_D,_reset_P)(a->ptr);
#endif
    Qs(Pst_args).a = 0;
    Qs(Pst_args).nc = -1;

    return 1;
}

static void
Qs(do_Ptrace)(QLA_D_Complex *r, int idx)
{
    Qx(QLA_D,_DiracPropagator) *a = &Qs(Pst_args).a[idx];
    int nc = Qs(Pst_args).nc;
    int is, ic;

    QLA_c_eq_r_plus_ir(*r, 0, 0);
    for (is = 0; is < QDP_Ns; is++) {
        for (ic = 0; ic < nc; ic++)  {
            QLA_c_peq_c(*r, Qx(QLA_D,_elem_P)(*a,ic,is,ic,is));
        }
    }
}

static int
Qs(q_P_trace)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    Qs(Pst_args).nc = QC(a);
#if QNc == 'N'
    Qs(Pst_args).a = Qx(QDP_D,_expose_P)(QC(a), a->ptr);
    QDP_D_C_eq_funci(r->ptr, Qs(do_Ptrace), *S->qss);
    Qx(QDP_D,_reset_P)(QC(a), a->ptr);
#else
    Qs(Pst_args).a = Qx(QDP_D,_expose_P)(a->ptr);
    QDP_D_C_eq_funci(r->ptr, Qs(do_Ptrace), *S->qss);
    Qx(QDP_D,_reset_P)(a->ptr);
#endif
    Qs(Pst_args).a = 0;
    Qs(Pst_args).nc = -1;

    return 1;
}

static int
Qs(q_P_set)(lua_State *L)
{
    Qs(mLatDirProp) *r = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 2, S, QC(r));

    CALL_QDP(L);
#if QNc == 'N'
    Qx(QDP_D,_P_eq_P)(QC(r), r->ptr, a->ptr, *S->qss);
#else
    Qx(QDP_D,_P_eq_P)(r->ptr, a->ptr, *S->qss);
#endif
    lua_pop(L, 1);

    return 1;
}

static int
Qs(q_P_add_P_)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirProp) *b = Qs(qlua_checkLatDirProp)(L, 2, S, QC(a));
    Qs(mLatDirProp) *c = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
#if QNc == 'N'
    Qx(QDP_D,_P_eq_P_plus_P)(QC(a), c->ptr, a->ptr, b->ptr, *S->qss);
#else
    Qx(QDP_D,_P_eq_P_plus_P)(c->ptr, a->ptr, b->ptr, *S->qss);
#endif

    return 1;
}

static int
Qs(q_P_sub_P_)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirProp) *b = Qs(qlua_checkLatDirProp)(L, 2, S, QC(a));
    Qs(mLatDirProp) *c = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
#if QNc == 'N'
    Qx(QDP_D,_P_eq_P_minus_P)(QC(a), c->ptr, a->ptr, b->ptr, *S->qss);
#else
    Qx(QDP_D,_P_eq_P_minus_P)(c->ptr, a->ptr, b->ptr, *S->qss);
#endif

    return 1;
}

static int
Qs(q_P_mul_r_)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_D_Real b = luaL_checknumber(L, 2);
    Qs(mLatDirProp) *c = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
#if QNc == 'N'
    Qx(QDP_D,_P_eq_r_times_P)(QC(a), c->ptr, &b, a->ptr, *S->qss);
#else
    Qx(QDP_D,_P_eq_r_times_P)(c->ptr, &b, a->ptr, *S->qss);
#endif

    return 1;
}

static int
Qs(q_r_mul_P_)(lua_State *L)
{
    QLA_D_Real a = luaL_checknumber(L, 1);
    Qs(mLatDirProp) *b = Qs(qlua_checkLatDirProp)(L, 2, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 2);
    Qs(mLatDirProp) *c = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(b));

    CALL_QDP(L);
#if QNc == 'N'
    Qx(QDP_D,_P_eq_r_times_P)(QC(b), c->ptr, &a, b->ptr, *S->qss);
#else
    Qx(QDP_D,_P_eq_r_times_P)(c->ptr, &a, b->ptr, *S->qss);
#endif

    return 1;
}

static int
Qs(q_P_mul_c_)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_D_Complex *b = qlua_checkComplex(L, 2);
    Qs(mLatDirProp) *c = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
#if QNc == 'N'
    Qx(QDP_D,_P_eq_c_times_P)(QC(a), c->ptr, b, a->ptr, *S->qss);
#else
    Qx(QDP_D,_P_eq_c_times_P)(c->ptr, b, a->ptr, *S->qss);
#endif

    return 1;
}

static int
Qs(q_c_mul_P_)(lua_State *L)
{
    QLA_D_Complex *a = qlua_checkComplex(L, 1);
    Qs(mLatDirProp) *b = Qs(qlua_checkLatDirProp)(L, 2, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 2);
    Qs(mLatDirProp) *c = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(b));

    CALL_QDP(L);
#if QNc == 'N'
    Qx(QDP_D,_P_eq_c_times_P)(QC(b), c->ptr, a, b->ptr, *S->qss);
#else
    Qx(QDP_D,_P_eq_c_times_P)(c->ptr, a, b->ptr, *S->qss);
#endif

    return 1;
}

static struct {
    int nc;
    QLA_D_Real *a;
    Qx(QLA_D,_DiracPropagator) *b;
} Qs(RPmul_args); /* YYY global state */

static void
Qs(do_RPmul)(Qx(QLA_D,_DiracPropagator) *r, int idx)
{
#if QNc == 'N'
    Qx(QLA_D,_P_eq_r_times_P)(Qs(RPmul_args).nc, r,
                              &Qs(RPmul_args).a[idx], &Qs(RPmul_args).b[idx]);
#else
    Qx(QLA_D,_P_eq_r_times_P)(r,
                              &Qs(RPmul_args).a[idx], &Qs(RPmul_args).b[idx]);
#endif
}

static void
Qs(X_P_eq_R_times_P)(Qx(QDP_D,_DiracPropagator) *r,
                     QDP_D_Real *a,
                     Qx(QDP_D,_DiracPropagator) *b,
                     QDP_Subset s,
                     int nc)
{
    Qs(RPmul_args).nc = nc;
    Qs(RPmul_args).a = QDP_D_expose_R(a);
#if QNc == 'N'
    Qs(RPmul_args).b = Qx(QDP_D,_expose_P)(nc, b);
    Qx(QDP_D,_P_eq_funci)(nc, r, Qs(do_RPmul), s);
    Qx(QDP_D,_reset_P)(nc, b);
#else
    Qs(RPmul_args).b = Qx(QDP_D,_expose_P)(b);
    Qx(QDP_D,_P_eq_funci)(r, Qs(do_RPmul), s);
    Qx(QDP_D,_reset_P)(b);
#endif
    QDP_D_reset_R(a);
    Qs(RPmul_args).a = 0;
    Qs(RPmul_args).b = 0;
    Qs(RPmul_args).nc = -1;
}

static int
Qs(q_R_mul_P_)(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirProp) *b = Qs(qlua_checkLatDirProp)(L, 2, S, -1);
    Qs(mLatDirProp) *r = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(b));

    CALL_QDP(L);
    Qs(X_P_eq_R_times_P)(r->ptr, a->ptr, b->ptr, *S->qss, QC(b));

    return 1;
}

static int
Qs(q_P_mul_R_)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, S);
    Qs(mLatDirProp) *r = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qs(X_P_eq_R_times_P)(r->ptr, b->ptr, a->ptr, *S->qss, QC(a));

    return 1;
}

static struct {
    int nc;
    QLA_D_Complex *a;
    Qx(QLA_D,_DiracPropagator) *b;
} Qs(CPmul_args); /* YYY global state */

static void
Qs(do_CPmul)(Qx(QLA_D,_DiracPropagator) *r, int idx)
{
#if QNc == 'N'
    Qx(QLA_D,_P_eq_c_times_P)(Qs(CPmul_args).nc, r,
                              &Qs(CPmul_args).a[idx], &Qs(CPmul_args).b[idx]);
#else
    Qx(QLA_D,_P_eq_c_times_P)(r,
                              &Qs(CPmul_args).a[idx], &Qs(CPmul_args).b[idx]);
#endif
}

static void
Qs(X_P_eq_C_times_P)(Qx(QDP_D,_DiracPropagator) *r,
                     QDP_D_Complex *a,
                     Qx(QDP_D,_DiracPropagator) *b,
                     QDP_Subset s,
                     int nc)
{
    Qs(CPmul_args).nc = nc;
    Qs(CPmul_args).a = QDP_D_expose_C(a);
#if QNc == 'N'
    Qs(CPmul_args).b = Qx(QDP_D,_expose_P)(nc, b);
    Qx(QDP_D,_P_eq_funci)(nc, r, Qs(do_CPmul), s);
    Qx(QDP_D,_reset_P)(nc, b);
#else
    Qs(CPmul_args).b = Qx(QDP_D,_expose_P)(b);
    Qx(QDP_D,_P_eq_funci)(r, Qs(do_CPmul), s);
    Qx(QDP_D,_reset_P)(b);
#endif
    QDP_D_reset_C(a);
    Qs(CPmul_args).a = 0;
    Qs(CPmul_args).b = 0;
    Qs(CPmul_args).nc = -1;
}

static int
Qs(q_C_mul_P_)(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirProp) *b = Qs(qlua_checkLatDirProp)(L, 2, S, -1);
    Qs(mLatDirProp) *r = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(b));

    CALL_QDP(L);
    Qs(X_P_eq_C_times_P)(r->ptr, a->ptr, b->ptr, *S->qss, QC(b));

    return 1;
}

static int
Qs(q_P_mul_C_)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2, S);
    Qs(mLatDirProp) *r = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);

    Qs(X_P_eq_C_times_P)(r->ptr, b->ptr, a->ptr, *S->qss, QC(a));

    return 1;
}

static int
Qs(q_P_neg)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirProp) *r = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
#if QNc == 'N'
    Qx(QDP_D,_P_meq_P)(QC(a), r->ptr, a->ptr, *S->qss);
#else
    Qx(QDP_D,_P_meq_P)(r->ptr, a->ptr, *S->qss);
#endif

    return 1;
}

static int
Qs(q_P_mul_P_)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirProp) *b = Qs(qlua_checkLatDirProp)(L, 2, S, QC(a));
    Qs(mLatDirProp) *c = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
#if QNc == 'N'
    Qx(QDP_D,_P_eq_P_times_P)(QC(a), c->ptr, a->ptr, b->ptr, *S->qss);
#else
    Qx(QDP_D,_P_eq_P_times_P)(c->ptr, a->ptr, b->ptr, *S->qss);
#endif

    return 1;
}

static int
Qs(q_P_mul_M_)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColMat) *b = Qs(qlua_checkLatColMat)(L, 2, S, QC(a));
    Qs(mLatDirProp) *c = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
#if QNc == 'N'
    Qx(QDP_D,_P_eq_P_times_M)(QC(a), c->ptr, a->ptr, b->ptr, *S->qss);
#else
    Qx(QDP_D,_P_eq_P_times_M)(c->ptr, a->ptr, b->ptr, *S->qss);
#endif

    return 1;
}

static int
Qs(q_M_mul_P_)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirProp) *b = Qs(qlua_checkLatDirProp)(L, 2, S, QC(a));
    Qs(mLatDirProp) *c = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
#if QNc == 'N'
    Qx(QDP_D,_P_eq_M_times_P)(QC(a), c->ptr, a->ptr, b->ptr, *S->qss);
#else
    Qx(QDP_D,_P_eq_M_times_P)(c->ptr, a->ptr, b->ptr, *S->qss);
#endif

    return 1;
}

static int
Qs(q_P_div_r_)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_D_Real b = 1 / luaL_checknumber(L, 2);
    Qs(mLatDirProp) *c = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
#if QNc == 'N'
    Qx(QDP_D,_P_eq_r_times_P)(QC(a), c->ptr, &b, a->ptr, *S->qss);
#else
    Qx(QDP_D,_P_eq_r_times_P)(c->ptr, &b, a->ptr, *S->qss);
#endif

    return 1;
}

static int
Qs(q_P_div_c_)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_D_Complex *b = qlua_checkComplex(L, 2);
    Qs(mLatDirProp) *c = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(a));
    double n = 1 / (QLA_real(*b) * QLA_real(*b) + QLA_imag(*b) * QLA_imag(*b));
    QLA_D_Complex s;

    CALL_QDP(L);
    QLA_real(s) = n * QLA_real(*b);
    QLA_imag(s) = -n * QLA_imag(*b);
#if QNc == 'N'
    Qx(QDP_D,_P_eq_c_times_P)(QC(a), c->ptr, &s, a->ptr, *S->qss);
#else
    Qx(QDP_D,_P_eq_c_times_P)(c->ptr, &s, a->ptr, *S->qss);
#endif

    return 1;
}

static int
Qs(q_P_dot_)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirProp) *b = Qs(qlua_checkLatDirProp)(L, 2, S, QC(a));
    mLatComplex *s = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
#if QNc == 'N'
    Qx(QDP_D,_C_eq_P_dot_P)(QC(a), s->ptr, a->ptr, b->ptr, *S->qss);
#else
    Qx(QDP_D,_C_eq_P_dot_P)(s->ptr, a->ptr, b->ptr, *S->qss);
#endif

    return 1;
}

static int
Qs(q_P_colors)(lua_State *L)
{
#if QNc == 'N'
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
#else
    Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
#endif
    lua_pushnumber(L, QC(a));

    return 1;
}

static int
Qs(q_P_copy)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirProp) *r = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
#if QNc == 'N'
    Qx(QDP_D,_P_eq_P)(QC(a), r->ptr, a->ptr, *S->qss);
#else
    Qx(QDP_D,_P_eq_P)(r->ptr, a->ptr, *S->qss);
#endif

    return 1;
}

static struct luaL_Reg Qs(mtLatDirProp)[] = {
    { "__tostring",        Qs(q_P_fmt)            },
    { "__gc",              Qs(q_P_gc)             },
    { "__index",           Qs(q_P_get)            },
    { "__newindex",        Qs(q_P_put)            },
    { "__unm",             Qs(q_P_neg)            },
    { "__add",             qlua_add               },
    { "__sub",             qlua_sub               },
    { "__mul",             qlua_mul               },
    { "__div",             qlua_div               },
    { "norm2",             Qs(q_P_norm2_)         },
    { "shift",             Qs(q_P_shift)          },
    { "conj",              Qs(q_P_conj)           },
    { "transpose",         Qs(q_P_trans)          },
    { "adjoin",            Qs(q_P_adjoin)         },
    { "spintranspose",     Qs(q_P_spintranspose)  },
    { "spintrace",         Qs(q_P_spintrace)      },
    { "trace",             Qs(q_P_trace)          },
    { "set",               Qs(q_P_set)            },
    { "colors",            Qs(q_P_colors)         },
    { "copy",              Qs(q_P_copy)           },
    /* "lattice" */
    /* "a-type" */
    { NULL,                NULL                   }
};

Qs(mLatDirProp) *
Qs(qlua_newLatDirProp)(lua_State *L, int Sidx, int nc)
{
#if QNc == 'N'
    Qx(QDP_D,_DiracPropagator) *v = Qx(QDP_D,_create_P)(nc);
#else
    Qx(QDP_D,_DiracPropagator) *v = Qx(QDP_D,_create_P)();
#endif
    Qs(mLatDirProp) *hdr;

    if (v == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
#if QNc == 'N'
        v = Qx(QDP_D,_create_P)(nc);
#else
        v = Qx(QDP_D,_create_P)();
#endif
        if (v == 0)
            luaL_error(L, "not enough memory (QDP_DiracPropagator" Qcolors ")");
    }
    hdr = lua_newuserdata(L, sizeof (Qs(mLatDirProp)));
    hdr->ptr = v;
#if QNc == 'N'
    hdr->nc = nc;
#endif
    qlua_createLatticeTable(L, Sidx, Qs(mtLatDirProp), Qs(qLatDirProp),
                            Qs(LatDirPropName));
    lua_setmetatable(L, -2);

    return hdr;
}

Qs(mLatDirProp) *
Qs(qlua_checkLatDirProp)(lua_State *L, int idx, mLattice *S, int nc)
{
    void *v = qlua_checkLatticeType(L, idx, Qs(qLatDirProp),
                                    Qs(LatDirPropName));
    Qs(mLatDirProp) *z;

    if (S) {
        mLattice *S1 = qlua_ObjLattice(L, idx);
        if (S1->id != S->id)
            luaL_error(L, "%s on a wrong lattice", Qs(LatDirPropName));
        lua_pop(L, 1);
    }
    z = (Qs(mLatDirProp) *)v;
#if QNc == 'N'
    if (nc != -1) {
        if (z->nc != nc)
            luaL_error(L, "Wrong number of colors");
    }
#endif

    return z;
}

static int
Qs(q_latdirprop_)(lua_State *L, mLattice *S, int nc, int off)
{
    switch (lua_gettop(L) - off) {
    case 1: {
        Qs(mLatDirProp) *v = Qs(qlua_newLatDirProp)(L, 1, nc);

        CALL_QDP(L);
#if QNc == 'N'
        Qx(QDP_D,_P_eq_zero)(nc, v->ptr, *S->qss);
#else
        Qx(QDP_D,_P_eq_zero)(v->ptr, *S->qss);
#endif
        
        return 1;
    }
    case 2: {
        Qs(mLatDirProp) *v = Qs(qlua_newLatDirProp)(L, 1, nc);
        
        switch (qlua_qtype(L, 2 + off)) {
        case Qs(qLatColMat): {
            Qs(mLatColMat) *w = Qs(qlua_checkLatColMat)(L, 2 + off, S, nc);
            mLatComplex *c = qlua_newLatComplex(L, 1);
            int ic, jc, ks;

            CALL_QDP(L);
#if QNc == 'N'
            Qx(QDP_D,_P_eq_zero)(nc, v->ptr, *S->qss);
#else
            Qx(QDP_D,_P_eq_zero)(v->ptr, *S->qss);
#endif
            for (ic = 0; ic < nc; ic++) {
                for (jc = 0; jc < nc; jc++) {
#if QNc == 'N'
                    Qx(QDP_D,_C_eq_elem_M)(nc, c->ptr, w->ptr, ic, jc, *S->qss);
                    for (ks = 0; ks < QDP_Ns; ks++)
                        Qx(QDP_D,_P_eq_elem_C)(nc, v->ptr, c->ptr,
                                               ic, ks, jc, ks, *S->qss);
#else
                    Qx(QDP_D,_C_eq_elem_M)(c->ptr, w->ptr, ic, jc, *S->qss);
                    for (ks = 0; ks < QDP_Ns; ks++)
                        Qx(QDP_D,_P_eq_elem_C)(v->ptr, c->ptr,
                                               ic, ks, jc, ks, *S->qss);
#endif
                }
            }
            lua_pop(L, 1);

            return 1;
        }
        default:
            break;
        }
    }
    case 3: {
        Qs(mLatDirFerm) *z = Qs(qlua_checkLatDirFerm)(L, 2 + off, S, nc);
        int d = qlua_checkdiracindex(L, 3);
        int c = qlua_checkcolorindex(L, 3, nc);
        Qs(mLatDirProp) *v = Qs(qlua_newLatDirProp)(L, 1, nc);

        CALL_QDP(L);
#if QNc == 'N'
        Qx(QDP_D,_P_eq_zero)(nc, v->ptr, *S->qss);
        Qx(QDP_D,_P_eq_diracvec_D)(nc, v->ptr, z->ptr, c, d, *S->qss);
#else
        Qx(QDP_D,_P_eq_zero)(v->ptr, *S->qss);
        Qx(QDP_D,_P_eq_diracvec_D)(v->ptr, z->ptr, c, d, *S->qss);
#endif
        return 1;
    }
    }
    return qlua_badconstr(L, "DiracPropagator" Qcolors);
}

#undef QNc
#undef Qcolors
#undef Qs
#undef Qx
#undef QC
