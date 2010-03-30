
static const char Qs(LatColVecName)[] = "lattice.ColorVector" Qcolors;

static int
Qs(q_V_fmt)(lua_State *L)
{
    char fmt[72];
    Qs(mLatColVec) *b = Qs(qlua_checkLatColVec)(L, 1, NULL, -1);

    sprintf(fmt, "QDP:ColorVector%d(%p)", QC(b), b->ptr);
    lua_pushstring(L, fmt);

    return 1;
}

static int
Qs(q_V_gc)(lua_State *L)
{
    Qs(mLatColVec) *b = Qs(qlua_checkLatColVec)(L, 1, NULL, -1);

    Qx(QDP_D,_destroy_V)(b->ptr);
    b->ptr = 0;

    return 0;
}

static int
Qs(q_V_get)(lua_State *L)
{
    switch (qlua_qtype(L, 2)) {
    case qTable: {
        Qs(mLatColVec) *V = Qs(qlua_checkLatColVec)(L, 1, NULL, -1);
        mLattice *S = qlua_ObjLattice(L, 1);
        int c = qlua_checkcolorindex(L, 2, QC(V));
        int *idx = qlua_latcoord(L, 2, S);
        int Sidx = lua_gettop(L);

        if (idx == NULL) {
            mLatComplex *r = qlua_newLatComplex(L, Sidx);

            CALL_QDP(L);
            Qx(QDP_D,_C_eq_elem_V)(r->ptr, V->ptr, c, *S->qss);
        } else {
            QLA_D_Complex *W = qlua_newComplex(L);
#if QNc == 'N'
            typedef QLA_DN_ColorVector(V->nc, Vtype);
#else
            typedef Qx(QLA_D,_ColorVector) Vtype;
#endif
            Vtype *locked;
            double zri[2];

            qlua_verifylatcoord(L, idx, S);
            CALL_QDP(L);
            locked = Qx(QDP_D,_expose_V)(V->ptr);

            if (QDP_node_number(idx) == QDP_this_node) {
                QLA_D_Complex *zz;

                zz = &Qx(QLA_D,_elem_V)(locked[QDP_index(idx)], c);
                zri[0] = QLA_real(*zz);
                zri[1] = QLA_imag(*zz);
            } else {
                zri[0] = 0;
                zri[1] = 0;
            }
            Qx(QDP_D,_reset_V)(V->ptr);
            QMP_sum_double_array(zri, 2);
            QLA_c_eq_r_plus_ir(*W, zri[0], zri[1]);
        }
        qlua_free(L, idx);
        return 1;
    }
    case qString:
        return  qlua_selflookup(L, 1, luaL_checkstring(L, 2));
    default:
        break;
    }
    return qlua_badindex(L, "ColorVector" Qcolors);
}

static int
Qs(q_V_put)(lua_State *L)
{
    Qs(mLatColVec) *V = Qs(qlua_checkLatColVec)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    int c = qlua_checkcolorindex(L, 2, QC(V));
    int *idx = qlua_latcoord(L, 2, S);

    if (idx == NULL) {
        mLatComplex *z = qlua_checkLatComplex(L, 3, S);
        
        CALL_QDP(L);
        Qx(QDP_D,_V_eq_elem_C)(V->ptr, z->ptr, c, *S->qss);
    } else {
        QLA_Complex *z = qlua_checkComplex(L, 3);
        qlua_verifylatcoord(L, idx, S);
#if QNc == 'N'
        typedef QLA_DN_ColorVector(V->nc, Vtype);
#else
        typedef Qx(QLA_D,_ColorVector) Vtype;
#endif
        Vtype *locked;
        QLA_Complex *zz;
        
        CALL_QDP(L);
        locked = Qx(QDP_D,_expose_V)(V->ptr);
        if (QDP_node_number(idx) == QDP_this_node) {
            zz = &Qx(QLA_D,_elem_V)(locked[QDP_index(idx)], c);
            QLA_c_eq_c(*zz, *z);
        }
        Qx(QDP_D,_reset_V)(V->ptr);
    }
    qlua_free(L, idx);
    return 0;
}

static int
Qs(q_V_dot_)(lua_State *L)
{
    Qs(mLatColVec) *a = Qs(qlua_checkLatColVec)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColVec) *b = Qs(qlua_checkLatColVec)(L, 2, S, QC(a));
    mLatComplex *s = qlua_newLatComplex(L, lua_gettop(L));

    CALL_QDP(L);
#if QNc == 'N'
    Qx(QDP_D,_C_eq_V_dot_V)(s->ptr, a->ptr, b->ptr, *S->qss);
#else
    Qx(QDP_D,_C_eq_V_dot_V)(s->ptr, a->ptr, b->ptr, *S->qss);
#endif

    return 1;
}

static int
Qs(q_V_add_V_)(lua_State *L)
{
    Qs(mLatColVec) *a = Qs(qlua_checkLatColVec)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColVec) *b = Qs(qlua_checkLatColVec)(L, 2, S, QC(a));
    Qs(mLatColVec) *c = Qs(qlua_newLatColVec)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_V_eq_V_plus_V)(c->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
Qs(q_V_sub_V_)(lua_State *L)
{
    Qs(mLatColVec) *a = Qs(qlua_checkLatColVec)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColVec) *b = Qs(qlua_checkLatColVec)(L, 2, S, QC(a));
    Qs(mLatColVec) *c = Qs(qlua_newLatColVec)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_V_eq_V_minus_V)(c->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
Qs(q_V_mul_r_)(lua_State *L)
{
    Qs(mLatColVec) *a = Qs(qlua_checkLatColVec)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_Real b = luaL_checknumber(L, 2);
    Qs(mLatColVec) *c = Qs(qlua_newLatColVec)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_V_eq_r_times_V)(c->ptr, &b, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_r_mul_V_)(lua_State *L)
{
    QLA_Real a = luaL_checknumber(L, 1);
    Qs(mLatColVec) *b = Qs(qlua_checkLatColVec)(L, 2, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 2);
    Qs(mLatColVec) *c = Qs(qlua_newLatColVec)(L, lua_gettop(L), QC(b));

    CALL_QDP(L);
    Qx(QDP_D,_V_eq_r_times_V)(c->ptr, &a, b->ptr, *S->qss);

    return 1;
}

static int
Qs(q_V_mul_c_)(lua_State *L)
{
    Qs(mLatColVec) *a = Qs(qlua_checkLatColVec)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    Qs(mLatColVec) *c = Qs(qlua_newLatColVec)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_V_eq_c_times_V)(c->ptr, b, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_c_mul_V_)(lua_State *L)
{
    QLA_Complex *a = qlua_checkComplex(L, 1);
    Qs(mLatColVec) *b = Qs(qlua_checkLatColVec)(L, 2, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 2);
    Qs(mLatColVec) *c = Qs(qlua_newLatColVec)(L, lua_gettop(L), QC(b));

    CALL_QDP(L);
    Qx(QDP_D,_V_eq_c_times_V)(c->ptr, a, b->ptr, *S->qss);

    return 1;
}

static struct {
    QLA_Real *a;
    int nc;
    void *b ; /* Qx(QLA_D,_ColorVector) */
} Qs(RVmul_args); /* YYY global state */

#if QNc == 'N'
static void
Qs(do_RVmul)(int nc, QLA_DN_ColorVector(nc, (*r)), int idx)
{
    QLA_DN_ColorVector(nc, (*b)) = Qs(RVmul_args).b;
    Qx(QLA_D,_V_eq_r_times_V)(nc, r, &Qs(RVmul_args).a[idx], &b[idx]);
}
#else
static void
Qs(do_RVmul)(Qx(QLA_D,_ColorVector) *r, int idx)
{
    Qx(QLA_D,_ColorVector) *b = Qs(RVmul_args).b;
    Qx(QLA_D,_V_eq_r_times_V)(r, &Qs(RVmul_args).a[idx], &b[idx]);
}
#endif

static void
Qs(X_V_eq_R_times_V)(int nc,
                     Qx(QDP_D,_ColorVector) *r,
                     QDP_Real *a,
                     Qx(QDP_D,_ColorVector) *b,
                     QDP_Subset s)
{
    Qs(RVmul_args).a = QDP_expose_R(a);
    Qs(RVmul_args).nc = nc;
    Qs(RVmul_args).b = Qx(QDP_D,_expose_V)(b);
    Qx(QDP_D,_V_eq_funci)(r, Qs(do_RVmul), s);
    Qx(QDP_D,_reset_V)(b);
    QDP_reset_R(a);
    Qs(RVmul_args).a = 0;
    Qs(RVmul_args).b = 0;
    Qs(RVmul_args).nc = -1;
}

static int
Qs(q_R_mul_V_)(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColVec) *b = Qs(qlua_checkLatColVec)(L, 2, S, -1);
    Qs(mLatColVec) *r = Qs(qlua_newLatColVec)(L, lua_gettop(L), QC(b));

    CALL_QDP(L);
    Qs(X_V_eq_R_times_V)(QC(b), r->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
Qs(q_V_mul_R_)(lua_State *L)
{
    Qs(mLatColVec) *a = Qs(qlua_checkLatColVec)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatReal *b = qlua_checkLatReal(L, 2, S);
    Qs(mLatColVec) *r = Qs(qlua_newLatColVec)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qs(X_V_eq_R_times_V)(QC(a), r->ptr, b->ptr, a->ptr, *S->qss);

    return 1;
}

static struct {
    int nc;
    QLA_Complex *a;
    void *b; /* Qx(QLA_D,_ColorVector)  */
} Qs(CVmul_args); /* YYY global state */

#if QNc == 'N'
static void
Qs(do_CVmul)(int nc, QLA_DN_ColorVector(nc, (*r)), int idx)
{
    QLA_DN_ColorVector(nc, (*b)) = Qs(CVmul_args).b;
    Qx(QLA_D,_V_eq_c_times_V)(nc, r, &Qs(CVmul_args).a[idx], &b[idx]);
}
#else
static void
Qs(do_CVmul)(Qx(QLA_D,_ColorVector) *r, int idx)
{
    Qx(QLA_D,_ColorVector) *b = Qs(CVmul_args).b;
    Qx(QLA_D,_V_eq_c_times_V)(r, &Qs(CVmul_args).a[idx], &b[idx]);
}
#endif

static void
Qs(X_V_eq_C_times_V)(int nc,
                     Qx(QDP_D,_ColorVector) *r,
                     QDP_Complex *a,
                     Qx(QDP_D,_ColorVector) *b,
                     QDP_Subset s)
{
    Qs(CVmul_args).nc = nc;
    Qs(CVmul_args).a = QDP_expose_C(a);
    Qs(CVmul_args).b = Qx(QDP_D,_expose_V)(b);
    Qx(QDP_D,_V_eq_funci)(r, Qs(do_CVmul), s);
    Qx(QDP_D,_reset_V)(b);
    QDP_reset_C(a);
    Qs(CVmul_args).a = 0;
    Qs(CVmul_args).b = 0;
    Qs(CVmul_args).nc = -1;
}

static int
Qs(q_C_mul_V_)(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColVec) *b = Qs(qlua_checkLatColVec)(L, 2, S, -1);
    Qs(mLatColVec) *r = Qs(qlua_newLatColVec)(L, lua_gettop(L), QC(b));

    CALL_QDP(L);
    Qs(X_V_eq_C_times_V)(QC(b), r->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
Qs(q_V_mul_C_)(lua_State *L)
{
    Qs(mLatColVec) *a = Qs(qlua_checkLatColVec)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatComplex *b = qlua_checkLatComplex(L, 2, S);
    Qs(mLatColVec) *r = Qs(qlua_newLatColVec)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qs(X_V_eq_C_times_V)(QC(a), r->ptr, b->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_V_norm2_)(lua_State *L)
{
    Qs(mLatColVec) *a = Qs(qlua_checkLatColVec)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_D_Real n;

    CALL_QDP(L);
    if (S->lss.mask) {
        Qs(mLatColVec) *b = Qs(qlua_newLatColVec)(L, lua_gettop(L), QC(a));
        Qx(QDP_D,_V_eq_V_mask_I)(b->ptr, a->ptr, S->lss.mask, *S->qss);
        Qx(QDP_D,_r_eq_norm2_V)(&n, b->ptr, *S->qss);
    } else {
        Qx(QDP_D,_r_eq_norm2_V)(&n, a->ptr, *S->qss);
    }
    lua_pushnumber(L, n);
    
    return 1;
}

static int
Qs(q_V_shift)(lua_State *L)
{
    Qs(mLatColVec) *a = Qs(qlua_checkLatColVec)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    QDP_Shift shift = qlua_checkShift(L, 2, S);
    QDP_ShiftDir dir = qlua_checkShiftDir(L, 3);
    Qs(mLatColVec) *r = Qs(qlua_newLatColVec)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_V_eq_sV)(r->ptr, a->ptr, shift, dir, *S->qss);

    return 1;
}

static int
Qs(q_V_conj)(lua_State *L)
{
    Qs(mLatColVec) *a = Qs(qlua_checkLatColVec)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColVec) *r = Qs(qlua_newLatColVec)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_V_eq_conj_V)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_V_set)(lua_State *L)
{
    Qs(mLatColVec) *r = Qs(qlua_checkLatColVec)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColVec) *a = Qs(qlua_checkLatColVec)(L, 2, S, QC(r));

    CALL_QDP(L);
    if (S->lss.mask)
        Qx(QDP_D,_V_eq_V_mask_I)(r->ptr, a->ptr, S->lss.mask, *S->qss);
    else
        Qx(QDP_D,_V_eq_V)(r->ptr, a->ptr, *S->qss);
    lua_pop(L, 2);

    return 1;
}

static int
Qs(q_V_neg)(lua_State *L)
{
    Qs(mLatColVec) *a = Qs(qlua_checkLatColVec)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColVec) *r = Qs(qlua_newLatColVec)(L, lua_gettop(L), QC(a));
    QLA_D_Real m1 = -1;

    CALL_QDP(L);
    Qx(QDP_D,_V_eq_r_times_V)(r->ptr, &m1, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_V_div_r_)(lua_State *L)
{
    Qs(mLatColVec) *a = Qs(qlua_checkLatColVec)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_D_Real b = 1 / luaL_checknumber(L, 2);
    Qs(mLatColVec) *c = Qs(qlua_newLatColVec)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_V_eq_r_times_V)(c->ptr, &b, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_V_div_c_)(lua_State *L)
{
    Qs(mLatColVec) *a = Qs(qlua_checkLatColVec)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_D_Complex *b = qlua_checkComplex(L, 2);
    Qs(mLatColVec) *c = Qs(qlua_newLatColVec)(L, lua_gettop(L), QC(a));
    double n = 1 / (QLA_real(*b) * QLA_real(*b) + QLA_imag(*b) * QLA_imag(*b));
    QLA_D_Complex s;

    CALL_QDP(L);
    QLA_real(s) = n * QLA_real(*b);
    QLA_imag(s) = -n * QLA_imag(*b);
    Qx(QDP_D,_V_eq_c_times_V)(c->ptr, &s, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_V_colors)(lua_State *L)
{
#if QNc == 'N'
    Qs(mLatColVec) *a = Qs(qlua_checkLatColVec)(L, 1, NULL, -1);
#else
    Qs(qlua_checkLatColVec)(L, 1, NULL, -1);
#endif
    lua_pushnumber(L, QC(a));

    return 1;
}

static int
Qs(q_V_copy)(lua_State *L)
{
    Qs(mLatColVec) *a = Qs(qlua_checkLatColVec)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColVec) *r = Qs(qlua_newLatColVec)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_V_eq_V)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static struct luaL_Reg Qs(mtLatColVec)[] = {
    { "__tostring",        Qs(q_V_fmt)    },
    { "__gc",              Qs(q_V_gc)     },
    { "__index",           Qs(q_V_get)    },
    { "__newindex",        Qs(q_V_put)    },
    { "__unm",             Qs(q_V_neg)    },
    { "__add",             qlua_add       },
    { "__sub",             qlua_sub       },
    { "__mul",             qlua_mul       },
    { "__div",             qlua_div       },
    { "norm2",             Qs(q_V_norm2_) },
    { "shift",             Qs(q_V_shift)  },
    { "conj",              Qs(q_V_conj)   },
    { "set",               Qs(q_V_set)    },
    { "colors",            Qs(q_V_colors) },
    { "copy",              Qs(q_V_copy)   },
    /* "lattice" */
    /* "a-type" */
    { NULL,                NULL           }
};

Qs(mLatColVec) *
Qs(qlua_newLatColVec)(lua_State *L, int Sidx, int nc)
{
#if QNc == 'N'
    Qx(QDP_D,_ColorVector) *v = Qx(QDP_D,_create_V)(nc);
#else
    Qx(QDP_D,_ColorVector) *v = Qx(QDP_D,_create_V)();
#endif
    Qs(mLatColVec) *hdr;

    if (v == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
#if QNc == 'N'
        v = Qx(QDP_D,_create_V)(nc);
#else
        v = Qx(QDP_D,_create_V)();
#endif
        if (v == 0)
            luaL_error(L, "not enough memory (QDP_ColorVector" Qcolors ")");
    }
    hdr = lua_newuserdata(L, sizeof (Qs(mLatColVec)));
    hdr->ptr = v;
#if QNc == 'N'
    hdr->nc = nc;
#endif
    qlua_createLatticeTable(L, Sidx, Qs(mtLatColVec), Qs(qLatColVec),
                            Qs(LatColVecName));
    lua_setmetatable(L, -2);

    return hdr;
}

Qs(mLatColVec) *
Qs(qlua_checkLatColVec)(lua_State *L, int idx, mLattice *S, int nc)
{
    void *v = qlua_checkLatticeType(L, idx, Qs(qLatColVec), Qs(LatColVecName));
    Qs(mLatColVec) *z;

    if (S) {
        mLattice *S1 = qlua_ObjLattice(L, idx);
        if (S1->id != S->id)
            luaL_error(L, "%s on a wrong lattice", Qs(LatColVecName));
        lua_pop(L, 1);
    }
    z = (Qs(mLatColVec) *)v;
#if QNc == 'N'
    if (nc != -1) {
        if (z->nc != nc)
            luaL_error(L, "Wrong number of colors");
    }
#endif

    return z;
}

static int
Qs(q_latcolvec_)(lua_State *L, mLattice *S, int nc, int off)
{
    switch (lua_gettop(L) - off) {
    case 1: {
        Qs(mLatColVec) *v = Qs(qlua_newLatColVec)(L, 1, nc);

        CALL_QDP(L);
        Qx(QDP_D,_V_eq_zero)(v->ptr, *S->qss);

        return 1;
    }
    case 2: {
        Qs(mLatColVec) *a = Qs(qlua_checkLatColVec)(L, 2 + off, S, nc);
        Qs(mLatColVec) *r = Qs(qlua_newLatColVec)(L, 1, nc);
        
        CALL_QDP(L);
        Qx(QDP_D,_V_eq_V)(r->ptr, a->ptr, *S->qss);
        
        return 1;
    }
    case 3: {
        mLatComplex *c = qlua_checkLatComplex(L, 2 + off, S);
        int a = qlua_checkcolorindex(L, 3 + off, nc);
        Qs(mLatColVec) *r = Qs(qlua_newLatColVec)(L, 1, nc);

        CALL_QDP(L);
        Qx(QDP_D,_V_eq_elem_C)(r->ptr, c->ptr, a, *S->qss);

        return 1;
    }
    }
    return qlua_badconstr(L, "ColorVector" Qcolors);
}

#undef QNc
#undef Qcolors
#undef Qs
#undef Qx
#undef QC
