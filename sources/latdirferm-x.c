
static const char Qs(LatDirFermName)[] = "lattice.DiracFermion" Qcolors;

static int
Qs(q_D_fmt)(lua_State *L)
{
    char fmt[72];
    Qs(mLatDirFerm) *b = Qs(qlua_checkLatDirFerm)(L, 1, NULL, -1);

    sprintf(fmt, "QDP:DiracFermion%d(%p)", QC(b), b->ptr);
    lua_pushstring(L, fmt);

    return 1;
}

static int
Qs(q_D_gc)(lua_State *L)
{
    Qs(mLatDirFerm) *b = Qs(qlua_checkLatDirFerm)(L, 1, NULL, -1);

    Qx(QDP_D,_destroy_D)(b->ptr);
    b->ptr = 0;

    return 0;
}

static int
Qs(q_D_get)(lua_State *L)
{
    switch (qlua_qtype(L, 2)) {
    case qTable: {
        Qs(mLatDirFerm) *V = Qs(qlua_checkLatDirFerm)(L, 1, NULL, -1);
        mLattice *S = qlua_ObjLattice(L, 1);
        int d = qlua_checkdiracindex(L, 2);
        int c = qlua_colorindex(L, 2, QC(V));
        int *idx = qlua_latcoord(L, 2, S);
        if (idx == NULL) {
            if (c == -1) {
                Qs(mLatColVec) *r = Qs(qlua_newLatColVec)(L, lua_gettop(L),
                                                          QC(V));

                CALL_QDP(L);
                Qx(QDP_D,_V_eq_colorvec_D)(r->ptr, V->ptr, d, *S->qss);
            } else {
                mLatComplex *r = qlua_newLatComplex(L, lua_gettop(L));

                CALL_QDP(L);
                Qx(QDP_D,_C_eq_elem_D)(r->ptr, V->ptr, c, d, *S->qss);
            }
        } else {
            if (c == -1) {
                qlua_free(L, idx);
                return qlua_badindex(L, "DiracFermion" Qcolors);
            } else {
                QLA_Complex *W = qlua_newComplex(L);
#if QNc == 'N'
                typedef QLA_DN_DiracFermion(V->nc, Vtype);
#else
                typedef Qx(QLA_D,_DiracFermion) Vtype;
#endif
                Vtype *locked;
                double zri[2];

                qlua_verifylatcoord(L, idx, S);
                CALL_QDP(L);
                locked = Qx(QDP_D,_expose_D)(V->ptr);
                if (QDP_node_number_L(S->lat, idx) == QDP_this_node) {
                    QLA_Complex *zz;
                    zz = &Qx(QLA_D,_elem_D)(locked[QDP_index_L(S->lat, idx)],
                                            c, d);
                    zri[0] = QLA_real(*zz);
                    zri[1] = QLA_imag(*zz);
                } else {
                    zri[0] = 0;
                    zri[1] = 0;
                }
                Qx(QDP_D,_reset_D)(V->ptr);
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
    return qlua_badindex(L, "DiracFermion" Qcolors);
}

static int
Qs(q_D_put)(lua_State *L)
{
    Qs(mLatDirFerm) *V = Qs(qlua_checkLatDirFerm)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    int d = qlua_checkdiracindex(L, 2);
    int c = qlua_colorindex(L, 2, QC(V));
    int *idx = qlua_latcoord(L, 2, S);

    if (idx == NULL) {
        if (c == -1) {
            Qs(mLatColVec) *z = Qs(qlua_checkLatColVec)(L, 3, S, QC(V));

            CALL_QDP(L);
            Qx(QDP_D,_D_eq_colorvec_V)(V->ptr, z->ptr, d, *S->qss);
        } else {
            mLatComplex *z = qlua_checkLatComplex(L, 3, S);

            CALL_QDP(L);
            Qx(QDP_D,_D_eq_elem_C)(V->ptr, z->ptr, c, d, *S->qss);
        }
    } else {
        if (c == -1) {
            qlua_free(L, idx);
            return qlua_badindex(L, "DiracFermion" Qcolors);
        } else {
            QLA_Complex *z = qlua_checkComplex(L, 3);
            qlua_verifylatcoord(L, idx, S);
            CALL_QDP(L);
#if QNc == 'N'
            typedef QLA_DN_DiracFermion(V->nc, Vtype);
#else
            typedef Qx(QLA_D,_DiracFermion) Vtype;
#endif
            Vtype *locked = Qx(QDP_D,_expose_D)(V->ptr);
            if (QDP_node_number_L(S->lat, idx) == QDP_this_node) {
                int ix = QDP_index_L(S->lat, idx);
                QLA_Complex *zz = &Qx(QLA_D,_elem_D)(locked[ix], c, d);

                QLA_c_eq_c(*zz, *z);
            }
            Qx(QDP_D,_reset_D)(V->ptr);
        }
    }
    qlua_free(L, idx);

    return 0;
}

static int
Qs(q_D_neg)(lua_State *L)
{
    Qs(mLatDirFerm) *a = Qs(qlua_checkLatDirFerm)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirFerm) *r = Qs(qlua_newLatDirFerm)(L, lua_gettop(L), QC(a));
    QLA_D_Real m1 = -1;

    CALL_QDP(L);
    Qx(QDP_D,_D_eq_r_times_D)(r->ptr, &m1, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_D_norm2_)(lua_State *L)
{
    Qs(mLatDirFerm) *a = Qs(qlua_checkLatDirFerm)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_D_Real n;

    CALL_QDP(L);
    if (S->lss.mask) {
        Qs(mLatDirFerm) *b = Qs(qlua_newLatDirFerm)(L, lua_gettop(L), QC(a));
        Qx(QDP_D,_D_eq_D_mask_I)(b->ptr, a->ptr, S->lss.mask, *S->qss);
        Qx(QDP_D,_r_eq_norm2_D)(&n, b->ptr, *S->qss);
    } else {
        Qx(QDP_D,_r_eq_norm2_D)(&n, a->ptr, *S->qss);
    }
    lua_pushnumber(L, n);
    
    return 1;
}

static int
Qs(q_D_shift)(lua_State *L)
{
    Qs(mLatDirFerm) *a = Qs(qlua_checkLatDirFerm)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    QDP_Shift shift = qlua_checkShift(L, 2, S);
    QDP_ShiftDir dir = qlua_checkShiftDir(L, 3);
    Qs(mLatDirFerm) *r = Qs(qlua_newLatDirFerm)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_D_eq_sD)(r->ptr, a->ptr, shift, dir, *S->qss);

    return 1;
}

static int
Qs(q_D_conj)(lua_State *L)
{
    Qs(mLatDirFerm) *a = Qs(qlua_checkLatDirFerm)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirFerm) *r = Qs(qlua_newLatDirFerm)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_D_eq_conj_D)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_D_gamma)(lua_State *L)
{
    Qs(mLatDirFerm) *f = Qs(qlua_checkLatDirFerm)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    int mu = qlua_gammaindex(L, 2);
    int n = qlua_gammabinary(L, 2);
    
    if (((n == -1) && (mu == -1)) || ((n != -1) && (mu != -1)))
        return qlua_badindex(L, "DiracFermion" Qcolors);
    if (n == -1) {
        Qs(mLatDirFerm) *r = Qs(qlua_newLatDirFerm)(L, lua_gettop(L), QC(f));

        CALL_QDP(L);
        if (mu < 5) {
            Qx(QDP_D,_D_eq_gamma_times_D)(r->ptr, f->ptr, 1 << mu, *S->qss);

        } else {
            Qx(QDP_D,_D_eq_gamma_times_D)(r->ptr, f->ptr, 15, *S->qss);
        }

        return 1;
    }
    if (mu == -1) {
        Qs(mLatDirFerm) *r = Qs(qlua_newLatDirFerm)(L, lua_gettop(L), QC(f));

        CALL_QDP(L);
        Qx(QDP_D,_D_eq_gamma_times_D)(r->ptr, f->ptr, n, *S->qss);
        
        return 1;
    }

    return luaL_error(L, "this could never happen");
}

static int
Qs(q_D_set)(lua_State *L)
{
    Qs(mLatDirFerm) *r = Qs(qlua_checkLatDirFerm)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirFerm) *a = Qs(qlua_checkLatDirFerm)(L, 2, S, QC(r));

    CALL_QDP(L);
    if (S->lss.mask)
        Qx(QDP_D,_D_eq_D_mask_I)(r->ptr, a->ptr, S->lss.mask, *S->qss);
    else
        Qx(QDP_D,_D_eq_D)(r->ptr, a->ptr, *S->qss);
    lua_pop(L, 2);

    return 0;
}

static int
Qs(q_D_dot_)(lua_State *L)
{
    Qs(mLatDirFerm) *a = Qs(qlua_checkLatDirFerm)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirFerm) *b = Qs(qlua_checkLatDirFerm)(L, 2, S, QC(a));
    mLatComplex *s = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
    Qx(QDP_D,_C_eq_D_dot_D)(s->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
Qs(q_D_add_D_)(lua_State *L)
{
    Qs(mLatDirFerm) *a = Qs(qlua_checkLatDirFerm)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirFerm) *b = Qs(qlua_checkLatDirFerm)(L, 2, S, QC(a));
    Qs(mLatDirFerm) *c = Qs(qlua_newLatDirFerm)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_D_eq_D_plus_D)(c->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
Qs(q_D_sub_D_)(lua_State *L)
{
    Qs(mLatDirFerm) *a = Qs(qlua_checkLatDirFerm)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirFerm) *b = Qs(qlua_checkLatDirFerm)(L, 2, S, QC(a));
    Qs(mLatDirFerm) *c = Qs(qlua_newLatDirFerm)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_D_eq_D_minus_D)(c->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
Qs(q_D_mul_r_)(lua_State *L)
{
    Qs(mLatDirFerm) *a = Qs(qlua_checkLatDirFerm)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_D_Real b = luaL_checknumber(L, 2);
    Qs(mLatDirFerm) *c = Qs(qlua_newLatDirFerm)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_D_eq_r_times_D)(c->ptr, &b, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_r_mul_D_)(lua_State *L)
{
    QLA_D_Real a = luaL_checknumber(L, 1);
    Qs(mLatDirFerm) *b = Qs(qlua_checkLatDirFerm)(L, 2, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 2);
    Qs(mLatDirFerm) *c = Qs(qlua_newLatDirFerm)(L, lua_gettop(L), QC(b));

    CALL_QDP(L);
    Qx(QDP_D,_D_eq_r_times_D)(c->ptr, &a, b->ptr, *S->qss);

    return 1;
}

static int
Qs(q_D_mul_c_)(lua_State *L)
{
    Qs(mLatDirFerm) *a = Qs(qlua_checkLatDirFerm)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_D_Complex *b = qlua_checkComplex(L, 2);
    Qs(mLatDirFerm) *c = Qs(qlua_newLatDirFerm)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_D_eq_c_times_D)(c->ptr, b, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_c_mul_D_)(lua_State *L)
{
    QLA_D_Complex *a = qlua_checkComplex(L, 1);
    Qs(mLatDirFerm) *b = Qs(qlua_checkLatDirFerm)(L, 2, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 2);
    Qs(mLatDirFerm) *c = Qs(qlua_newLatDirFerm)(L, lua_gettop(L), QC(b));

    CALL_QDP(L);
    Qx(QDP_D,_D_eq_c_times_D)(c->ptr, a, b->ptr, *S->qss);

    return 1;
}

static struct {
    int nc;
    QLA_D_Real *a;
    void *b; /* Qx(QLA_D,_DiracFermion) *b; */
} Qs(RDmul_args); /* YYY global state */


#if QNc == 'N'
static void
Qs(do_RDmul)(int nc, QLA_DN_DiracFermion(nc, (*r)), int idx)
{
    QLA_DN_DiracFermion(nc, (*b)) = Qs(RDmul_args).b;
    Qx(QLA_D,_D_eq_r_times_D)(nc, r, &Qs(RDmul_args).a[idx], &b[idx]);
}
#else
static void
Qs(do_RDmul)(Qx(QLA_D,_DiracFermion) *r, int idx)
{
    Qx(QLA_D,_DiracFermion) *b = Qs(RDmul_args).b;
    Qx(QLA_D,_D_eq_r_times_D)(r, &Qs(RDmul_args).a[idx], &b[idx]);
}
#endif

static void
Qs(X_D_eq_R_times_D)(int nc,
                     Qx(QDP_D,_DiracFermion) *r,
                     QDP_D_Real *a,
                     Qx(QDP_D,_DiracFermion) *b,
                     QDP_Subset s)
{
    Qs(RDmul_args).nc = nc;
    Qs(RDmul_args).a = QDP_D_expose_R(a);
    Qs(RDmul_args).b = Qx(QDP_D,_expose_D)(b);
    Qx(QDP_D,_D_eq_funci)(r, Qs(do_RDmul), s);
    Qx(QDP_D,_reset_D)(b);
    QDP_D_reset_R(a);
    Qs(RDmul_args).a = 0;
    Qs(RDmul_args).b = 0;
    Qs(RDmul_args).nc = -1;
}

static int
Qs(q_R_mul_D_)(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirFerm) *b = Qs(qlua_checkLatDirFerm)(L, 2, S, -1);
    Qs(mLatDirFerm) *c = Qs(qlua_newLatDirFerm)(L, lua_gettop(L), QC(b));

    CALL_QDP(L);
    Qs(X_D_eq_R_times_D)(QC(b), c->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
Qs(q_D_mul_R_)(lua_State *L)
{
    Qs(mLatDirFerm) *a = Qs(qlua_checkLatDirFerm)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, S);
    Qs(mLatDirFerm) *c = Qs(qlua_newLatDirFerm)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qs(X_D_eq_R_times_D)(QC(a), c->ptr, b->ptr, a->ptr, *S->qss);

    return 1;
}

static struct {
    int nc;
    QLA_D_Complex *a;
    void *b; /* Qx(QLA_D,_DiracFermion) *b; */
} Qs(CDmul_args); /* YYY global state */

#if QNc == 'N'
static void
Qs(do_CDmul)(int nc, QLA_DN_DiracFermion(nc, (*r)), int idx)
{
    QLA_DN_DiracFermion(nc, (*b)) = Qs(CDmul_args).b;
    Qx(QLA_D,_D_eq_c_times_D)(nc, r, &Qs(CDmul_args).a[idx], &b[idx]);
}
#else
static void
Qs(do_CDmul)(Qx(QLA_D,_DiracFermion) *r, int idx)
{
    Qx(QLA_D,_DiracFermion) *b = Qs(CDmul_args).b;
    Qx(QLA_D,_D_eq_c_times_D)(r, &Qs(CDmul_args).a[idx], &b[idx]);
}
#endif

static void
Qs(X_D_eq_C_times_D)(int nc,
                     Qx(QDP_D,_DiracFermion) *r,
                     QDP_D_Complex *a,
                     Qx(QDP_D,_DiracFermion) *b,
                     QDP_Subset s)
{
    Qs(CDmul_args).nc = nc;
    Qs(CDmul_args).a = QDP_D_expose_C(a);
    Qs(CDmul_args).b = Qx(QDP_D,_expose_D)(b);
    Qx(QDP_D,_D_eq_funci)(r, Qs(do_CDmul), s);
    Qx(QDP_D,_reset_D)(b);
    QDP_D_reset_C(a);
    Qs(CDmul_args).a = 0;
    Qs(CDmul_args).b = 0;
    Qs(CDmul_args).nc = -1;
}

static int
Qs(q_C_mul_D_)(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirFerm) *b = Qs(qlua_checkLatDirFerm)(L, 2, S, -1);
    Qs(mLatDirFerm) *c = Qs(qlua_newLatDirFerm)(L, lua_gettop(L), QC(b));

    CALL_QDP(L);
    Qs(X_D_eq_C_times_D)(QC(b), c->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
Qs(q_D_mul_C_)(lua_State *L)
{
    Qs(mLatDirFerm) *a = Qs(qlua_checkLatDirFerm)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2, S);
    Qs(mLatDirFerm) *c = Qs(qlua_newLatDirFerm)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qs(X_D_eq_C_times_D)(QC(a), c->ptr, b->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_M_mul_D_)(lua_State *L)
{
    Qs(mLatColMat) *a = Qs(qlua_checkLatColMat)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirFerm) *b = Qs(qlua_checkLatDirFerm)(L, 2, S, QC(a));
    Qs(mLatDirFerm) *c = Qs(qlua_newLatDirFerm)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_D_eq_M_times_D)(c->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
Qs(q_D_div_r_)(lua_State *L)
{
    Qs(mLatDirFerm) *a = Qs(qlua_checkLatDirFerm)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_D_Real b = 1 / luaL_checknumber(L, 2);
    Qs(mLatDirFerm) *c = Qs(qlua_newLatDirFerm)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_D_eq_r_times_D)(c->ptr, &b, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_D_div_c_)(lua_State *L)
{
    Qs(mLatDirFerm) *a = Qs(qlua_checkLatDirFerm)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_D_Complex *b = qlua_checkComplex(L, 2);
    Qs(mLatDirFerm) *c = Qs(qlua_newLatDirFerm)(L, lua_gettop(L), QC(a));
    double n = 1 / (QLA_real(*b) * QLA_real(*b) + QLA_imag(*b) * QLA_imag(*b));
    QLA_D_Complex s;

    CALL_QDP(L);
    QLA_real(s) = n * QLA_real(*b);
    QLA_imag(s) = -n * QLA_imag(*b);
    Qx(QDP_D,_D_eq_c_times_D)(c->ptr, &s, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_D_colors)(lua_State *L)
{
#if QNc == 'N'
    Qs(mLatDirFerm) *a = Qs(qlua_checkLatDirFerm)(L, 1, NULL, -1);
#else
    Qs(qlua_checkLatDirFerm)(L, 1, NULL, -1);
#endif
    lua_pushnumber(L, QC(a));

    return 1;
}

static int
Qs(q_D_copy)(lua_State *L)
{
    Qs(mLatDirFerm) *a = Qs(qlua_checkLatDirFerm)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirFerm) *r = Qs(qlua_newLatDirFerm)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_D_eq_D)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static struct luaL_Reg Qs(mtLatDirFerm)[] = {
    { "__tostring",        Qs(q_D_fmt)     },
    { "__gc",              Qs(q_D_gc)      },
    { "__index",           Qs(q_D_get)     },
    { "__newindex",        Qs(q_D_put)     },
    { "__unm",             Qs(q_D_neg)     },
    { "__add",             qlua_add        },
    { "__sub",             qlua_sub        },
    { "__mul",             qlua_mul        },
    { "__div",             qlua_div        },
    { "norm2",             Qs(q_D_norm2_)  },
    { "shift",             Qs(q_D_shift)   },
    { "conj",              Qs(q_D_conj)    },
    { "gamma",             Qs(q_D_gamma)   },
    { "set",               Qs(q_D_set)     },
    { "colors",            Qs(q_D_colors)  },
    { "copy",              Qs(q_D_copy)    },
    /* "lattice" */
    /* "a-type" */
    { NULL,                NULL }
};

Qs(mLatDirFerm) *
Qs(qlua_newLatDirFerm)(lua_State *L, int Sidx, int nc)
{
    mLattice *S = qlua_checkLattice(L, Sidx);
#if QNc == 'N'
    Qx(QDP_D,_DiracFermion) *v = Qx(QDP_D,_create_D_L)(nc, S->lat);
#else
    Qx(QDP_D,_DiracFermion) *v = Qx(QDP_D,_create_D_L)(S->lat);
#endif
    Qs(mLatDirFerm) *hdr;

    if (v == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
#if QNc == 'N'
        v = Qx(QDP_D,_create_D_L)(nc, S->lat);
#else
        v = Qx(QDP_D,_create_D_L)(S->lat);
#endif
        if (v == 0)
            luaL_error(L, "not enough memory (QDP_DiracFermion" Qcolors ")");
    }
    hdr = lua_newuserdata(L, sizeof (Qs(mLatDirFerm)));
    hdr->ptr = v;
#if QNc == 'N'
    hdr->nc = nc;
#endif
    qlua_createLatticeTable(L, Sidx, Qs(mtLatDirFerm), Qs(qLatDirFerm),
                            Qs(LatDirFermName));
    lua_setmetatable(L, -2);

    return hdr;
}

Qs(mLatDirFerm) *
Qs(qlua_checkLatDirFerm)(lua_State *L, int idx, mLattice *S, int nc)
{
    void *v = qlua_checkLatticeType(L, idx, Qs(qLatDirFerm),
                                    Qs(LatDirFermName));
    Qs(mLatDirFerm) *z;

    if (S) {
        mLattice *S1 = qlua_ObjLattice(L, idx);
        if (S1->id != S->id)
            luaL_error(L, "%s on a wrong lattice", Qs(LatDirFermName));
        lua_pop(L, 1);
    }
    z = (Qs(mLatDirFerm) *)v;
#if QNc == 'N'
    if (nc != -1) {
        if (z->nc != nc)
            luaL_error(L, "Wrong number of colors");
    }
#endif

    return z;
}

static int
Qs(q_latdirferm_)(lua_State *L, mLattice *S, int nc, int off)
{
    switch (lua_gettop(L) - off) {
    case 1: {
        Qs(mLatDirFerm) *v = Qs(qlua_newLatDirFerm)(L, 1, nc);

        CALL_QDP(L);
        Qx(QDP_D,_D_eq_zero)(v->ptr, *S->qss);

        return 1;
    }
    case 3: {
        switch (qlua_qtype(L, 2 + off)) {
        case qLatComplex: {
            mLatComplex *z = qlua_checkLatComplex(L, 2 + off, S);
            int c = qlua_checkcolorindex(L, 3 + off, nc);
            int d = qlua_checkdiracindex(L, 3 + off);
            Qs(mLatDirFerm) *v = Qs(qlua_newLatDirFerm)(L, 1, nc);

            CALL_QDP(L);
            Qx(QDP_D,_D_eq_zero)(v->ptr, *S->qss);
            Qx(QDP_D,_D_eq_elem_C)(v->ptr, z->ptr, c, d, *S->qss);

            return 1;
        }
        case Qs(qLatColVec): {
            Qs(mLatColVec) *w = Qs(qlua_checkLatColVec)(L, 2 + off, S, nc);
            int d = qlua_checkdiracindex(L, 3 + off);
            Qs(mLatDirFerm) *v = Qs(qlua_newLatDirFerm)(L, 1, nc);

            CALL_QDP(L);
            Qx(QDP_D,_D_eq_zero)(v->ptr, *S->qss);
            Qx(QDP_D,_D_eq_colorvec_V)(v->ptr, w->ptr, d, *S->qss);
            return 1;
        }
        default:
            break;
        }
        break;
    }
    }
    return qlua_badconstr(L, "DiracFermion" Qcolors);
}

static const QLUA_Op2 Qs(ops)[] = {
    { qlua_add_table, Qs(qLatDirFerm),  Qs(qLatDirFerm),  Qs(q_D_add_D_) },
    { qlua_sub_table, Qs(qLatDirFerm),  Qs(qLatDirFerm),  Qs(q_D_sub_D_) },
    { qlua_mul_table, qReal,            Qs(qLatDirFerm),  Qs(q_r_mul_D_) },
    { qlua_mul_table, Qs(qLatDirFerm),  qReal,            Qs(q_D_mul_r_) },
    { qlua_mul_table, qComplex,         Qs(qLatDirFerm),  Qs(q_c_mul_D_) },
    { qlua_mul_table, Qs(qLatDirFerm),  qComplex,         Qs(q_D_mul_c_) },
    { qlua_mul_table, qLatReal,         Qs(qLatDirFerm),  Qs(q_R_mul_D_) },
    { qlua_mul_table, Qs(qLatDirFerm),  qLatReal,         Qs(q_D_mul_R_) },
    { qlua_mul_table, qLatComplex,      Qs(qLatDirFerm),  Qs(q_C_mul_D_) },
    { qlua_mul_table, Qs(qLatDirFerm),  qLatComplex,      Qs(q_D_mul_C_) },
    { qlua_mul_table, Qs(qLatColMat),   Qs(qLatDirFerm),  Qs(q_M_mul_D_) },
    { qlua_div_table, Qs(qLatDirFerm),  qReal,            Qs(q_D_div_r_) },
    { qlua_div_table, Qs(qLatDirFerm),  qComplex,         Qs(q_D_div_c_) },
    { NULL,           qNoType,          qNoType,          NULL           }
};

#undef QNc
#undef Qcolors
#undef Qs
#undef Qx
#undef QC
