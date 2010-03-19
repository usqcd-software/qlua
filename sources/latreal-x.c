const char Qa2(LatReal,Name)[] = "lattice.Real";

static int
Qa2(q_R,_fmt)(lua_State *L)
{
    char fmt[72];
    Qa1(mLatReal) *b = Qa1(qlua_checkLatReal)(L, 1, NULL);

    sprintf(fmt, "LatReal%c(%p)", Qx, b->ptr);
    lua_pushstring(L, fmt);

    return 1;
}

static int
Qa2(q_R,_gc)(lua_State *L)
{
    Qa1(mLatReal) *b = Qa1(qlua_checkLatReal)(L, 1, NULL);

    Qa2(QDP_,_destroy_R)(b->ptr);
    b->ptr = 0;

    return 0;
}

static int
Qa2(q_R,_neg)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa1(mLatReal) *res = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_R_eqm_R)(res->ptr, a->ptr, *S->qss);
    return 1;
}

static int
Qa2(q_R,_sum)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    int argc = lua_gettop(L);
    mLattice *S = qlua_ObjLattice(L, 1);
    
    switch (argc) {
    case 1: {
        Qa2(QLA_,_Real) sum;

        CALL_QDP(L);
        Qa2(QDP_,_r_eq_sum_R)(&sum, a->ptr, *S->qss);
        lua_pushnumber(L, sum);

        return 1;
    }
    case 2: {
        mLatMulti *m = qlua_checkLatMulti(L, 2, S);
        int size = m->size;
        QLA_Int *ii = m->idx;
        mVecReal *r = qlua_newVecReal(L, size);
        int k;
        Qa2(QLA_,_Real) *xx;
        
        for (k = 0; k < size; k++)
            r->val[k] = 0;

        CALL_QDP(L);
        xx = Qa2(QDP_,_expose_R)(a->ptr);
        for (k = 0; k < QDP_sites_on_node; k++, xx++, ii++) {
            int t = *ii;
            if ((t < 0) || (t >= size))
                continue;
            r->val[t] += *xx;
        }
        Qa2(QDP_,_reset_R)(a->ptr);
        QMP_sum_double_array(r->val, size);

        return 1;
    }
    }
    return luaL_error(L, "bad arguments for Real:sum()");
}

static int
Qa2(q_R,_norm2)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa2(QLA_,_Real) n;

    CALL_QDP(L);
    Qa2(QDP_,_r_eq_norm2_R)(&n, a->ptr, *S->qss);
    lua_pushnumber(L, n);

    return 1;
}

static int
Qa2(q_R,_shift)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    QDP_Shift shift = qlua_checkShift(L, 2, S);
    QDP_ShiftDir dir = qlua_checkShiftDir(L, 3);
    Qa1(mLatReal) *r = Qa1(qlua_newLatReal)(L, Sidx);

    CALL_QDP(L);
    Qa2(QDP_,_R_eq_sR)(r->ptr, a->ptr, shift, dir, *S->qss);

    return 1;
}

static int
Qa2(q_R,_sin)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa1(mLatReal) *r = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_R_eq_sin_R)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qa2(q_R,_cos)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa1(mLatReal) *r = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_R_eq_cos_R)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qa2(q_R,_tan)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa1(mLatReal) *r = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_R_eq_tan_R)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qa2(q_R,_asin)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa1(mLatReal) *r = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_R_eq_asin_R)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qa2(q_R,_acos)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa1(mLatReal) *r = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_R_eq_acos_R)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qa2(q_R,_atan)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa1(mLatReal) *r = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_R_eq_atan_R)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qa2(q_R,_sqrt)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa1(mLatReal) *r = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_R_eq_sqrt_R)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qa2(q_R,_abs)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa1(mLatReal) *r = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_R_eq_fabs_R)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qa2(q_R,_exp)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa1(mLatReal) *r = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_R_eq_exp_R)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qa2(q_R,_expi)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa1(mLatComplex) *r = Qa1(qlua_newLatComplex)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_C_eq_cexpi_R)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qa2(q_R,_log)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa1(mLatReal) *r = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_R_eq_log_R)(r->ptr, a->ptr, *S->qss);

    return 1;

}

static int
Qa2(q_R,_sign)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa1(mLatReal) *r = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_R_eq_sign_R)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qa2(q_R,_floor)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa1(mLatReal) *r = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_R_eq_floor_R)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qa2(q_R,_ceil)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa1(mLatReal) *r = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_R_eq_ceil_R)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qa2(q_R,_sinh)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa1(mLatReal) *r = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_R_eq_sinh_R)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qa2(q_R,_cosh)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa1(mLatReal) *r = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_R_eq_cosh_R)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qa2(q_R,_tanh)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa1(mLatReal) *r = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_R_eq_tanh_R)(r->ptr, a->ptr, *S->qss);

    return 1;

}

static int
Qa2(q_R,_log10)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa1(mLatReal) *r = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_R_eq_log10_R)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qa2(q_R,_trunc)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_I_eq_trunc_R)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qa2(q_R,_round)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatInt *r = qlua_newLatInt(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_I_eq_round_R)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qa2(q_R,_set)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa1(mLatReal) *r = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_R_eq_R)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qa2(q_R,_get)(lua_State *L)
{
    switch (qlua_qtype(L, 2)) {
    case qTable: {
        Qa1(mLatReal) *V = Qa1(qlua_checkLatReal)(L, 1, NULL);
        mLattice *S = qlua_ObjLattice(L, 1);
        Qa2(QLA_,_Real) *locked;
        int *idx = 0;
        double z;

        idx = qlua_checklatcoord(L, 2, S);
        CALL_QDP(L);
        locked = Qa2(QDP_,_expose_R)(V->ptr);
        if (QDP_node_number(idx) == QDP_this_node) {
            z = Qa2(QLA_,_elem_R)(locked[QDP_index(idx)]);
        } else {
            z = 0;
        }
        Qa2(QDP_,_reset_R)(V->ptr);
        qlua_free(L, idx);
        QMP_sum_double(&z);
        lua_pushnumber(L, z);

        return 1;
    }
    case qString:
        return qlua_selflookup(L, 1, luaL_checkstring(L, 2));
    default:
        break;
    }
    return qlua_badindex(L, "Real");
}


static int
Qa2(q_R,_put)(lua_State *L)
{
    Qa1(mLatReal) *V = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa2(QLA_,_Real) *locked;
    int *idx = 0;
    double z = luaL_checknumber(L, 3);

    idx = qlua_checklatcoord(L, 2, S);
    CALL_QDP(L);
    locked = Qa2(QDP_,_expose_R)(V->ptr);
    if (QDP_node_number(idx) == QDP_this_node) {
        Qa2(QLA_,_elem_R)(locked[QDP_index(idx)]) = z;
    }
    Qa2(QDP_,_reset_R)(V->ptr);
    qlua_free(L, idx);

    return 0;
}

static int
Qab(q_R,_add_R)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa1(mLatReal) *b = Qa1(qlua_checkLatReal)(L, 2, S);
    Qa1(mLatReal) *c = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_R_eq_R_plus_R)(c->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static struct {
    Qa2(QLA_,_Real) a;
    Qa2(QLA_,_Real) *b;
} Qa1(Rop_args); /* YYY global state */

static void
Qa1(do_rRadd)(Qa2(QLA_,_Real) *r, int idx)
{
    *r = Qa1(Rop_args).a + Qa1(Rop_args).b[idx];
}


static int
Qa1(q_r_add_R)(lua_State *L)
{
    double a = luaL_checknumber(L, 1);
    Qa1(mLatReal) *b = Qa1(qlua_checkLatReal)(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    Qa1(mLatReal) *c = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa1(Rop_args).a = a;
    Qa1(Rop_args).b = Qa2(QDP_,_expose_R)(b->ptr);
    Qa2(QDP_,_R_eq_funci)(c->ptr, Qa1(do_rRadd), *S->qss);
    Qa2(QDP_,_reset_R)(b->ptr);
    Qa1(Rop_args).a = 0;
    Qa1(Rop_args).b = 0;

    return 1;
}

static int
Qa2(q_R,_add_r)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    double b = luaL_checknumber(L, 2);
    Qa1(mLatReal) *c = Qa1(qlua_newLatReal)(L, lua_gettop(L));


    CALL_QDP(L);
    Qa1(Rop_args).a = b;
    Qa1(Rop_args).b = Qa2(QDP_,_expose_R)(a->ptr);
    Qa2(QDP_,_R_eq_funci)(c->ptr, Qa1(do_rRadd), *S->qss);
    Qa2(QDP_,_reset_R)(a->ptr);
    Qa1(Rop_args).a = 0;
    Qa1(Rop_args).b = 0;

    return 1;
}

static int
Qab(q_R,_sub_R)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa1(mLatReal) *b = Qa1(qlua_checkLatReal)(L, 2, S);
    Qa1(mLatReal) *c = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_R_eq_R_minus_R)(c->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static void
Qa1(do_rRsub)(Qa2(QLA_,_Real) *r, int idx)
{
    *r = Qa1(Rop_args).a - Qa1(Rop_args).b[idx];
}

static int
Qa1(q_r_sub_R)(lua_State *L)
{
    double a = luaL_checknumber(L, 1);
    Qa1(mLatReal) *b = Qa1(qlua_checkLatReal)(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    Qa1(mLatReal) *c = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa1(Rop_args).a = a;
    Qa1(Rop_args).b = Qa2(QDP_,_expose_R)(b->ptr);
    Qa2(QDP_,_R_eq_funci)(c->ptr, Qa1(do_rRsub), *S->qss);
    Qa2(QDP_,_reset_R)(b->ptr);
    Qa1(Rop_args).a = 0;
    Qa1(Rop_args).b = 0;

    return 1;
}

static int
Qa2(q_R,_sub_r)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    double b = luaL_checknumber(L, 2);
    Qa1(mLatReal) *c = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa1(Rop_args).a = -b;
    Qa1(Rop_args).b = Qa2(QDP_,_expose_R)(a->ptr);
    Qa2(QDP_,_R_eq_funci)(c->ptr, Qa1(do_rRadd), *S->qss);
    Qa2(QDP_,_reset_R)(a->ptr);
    Qa1(Rop_args).a = 0;
    Qa1(Rop_args).b = 0;

    return 1;

}

static int
Qab(q_R,_mul_R)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa1(mLatReal) *b = Qa1(qlua_checkLatReal)(L, 2, S);
    Qa1(mLatReal) *c = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_R_eq_R_times_R)(c->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
Qa1(q_r_mul_R)(lua_State *L)
{
    Qa2(QLA_,_Real) b = luaL_checknumber(L, 1);
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    Qa1(mLatReal) *c = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_R_eq_r_times_R)(c->ptr, &b, a->ptr, *S->qss);

    return 1;
}

static int
Qa2(q_R,_mul_r)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa2(QLA_,_Real) b = luaL_checknumber(L, 2);
    Qa1(mLatReal) *c = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_R_eq_r_times_R)(c->ptr, &b, a->ptr, *S->qss);

    return 1;
}

static int
Qab(q_R,_div_R)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa1(mLatReal) *b = Qa1(qlua_checkLatReal)(L, 2, S);
    Qa1(mLatReal) *c = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_R_eq_R_divide_R)(c->ptr, a->ptr, b->ptr, *S->qss);

    return 1;
}

static int
Qa2(q_R,_div_r)(lua_State *L)
{
    Qa1(mLatReal) *a = Qa1(qlua_checkLatReal)(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa2(QLA_,_Real) b = 1/luaL_checknumber(L, 2);
    Qa1(mLatReal) *c = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_R_eq_r_times_R)(c->ptr, &b, a->ptr, *S->qss);

    return 1;
}

static void
Qa1(do_rRdiv)(Qa2(QLA_,_Real) *r, int idx)
{
    *r = Qa1(Rop_args).a / Qa1(Rop_args).b[idx];
}

static int
Qa1(q_r_div_R)(lua_State *L)
{
    double a = luaL_checknumber(L, 1);
    Qa1(mLatReal) *b = Qa1(qlua_checkLatReal)(L, 2, NULL);
    mLattice *S = qlua_ObjLattice(L, 2);
    Qa1(mLatReal) *c = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);

    Qa1(Rop_args).a = a;
    Qa1(Rop_args).b = Qa2(QDP_,_expose_R)(b->ptr);
    Qa2(QDP_,_R_eq_funci)(c->ptr, Qa1(do_rRdiv), *S->qss);
    Qa2(QDP_,_reset_R)(b->ptr);
    Qa1(Rop_args).a = 0;
    Qa1(Rop_args).b = 0;

    return 1;
}

static int
Qa1(q_latreal)(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);

    switch (lua_gettop(L)) {
    case 1: {
        Qa1(mLatReal) *v = Qa1(qlua_newLatReal)(L, 1);

        CALL_QDP(L);
        Qa2(QDP_,_R_eq_zero)(v->ptr, *S->qss);

        return 1;
    }
    case 2:
        switch (qlua_qtype(L, 2)) {
        case qReal: {
            Qa2(QLA_,_Real) d = luaL_checknumber(L, 2);
            Qa1(mLatReal) *v = Qa1(qlua_newLatReal)(L, 1);
            
            CALL_QDP(L);
            Qa2(QDP_,_R_eq_r)(v->ptr, &d, *S->qss);
            
            return 1;
        }
        case qLatInt: {
            mLatInt *d = qlua_checkLatInt(L, 2, S);
            Qa1(mLatReal) *v = Qa1(qlua_newLatReal)(L, 1);
            
            CALL_QDP(L);
            Qa2(QDP_,_R_eq_I)(v->ptr, d->ptr, *S->qss);
            
            return 1;
        }
        case Qa1(qLatReal): { /* no mixing constructor */
            Qa1(mLatReal) *d = Qa1(qlua_checkLatReal)(L, 2, S);
            Qa1(mLatReal) *v = Qa1(qlua_newLatReal)(L, 1);
            
            CALL_QDP(L);
            Qa2(QDP_,_R_eq_R)(v->ptr, d->ptr, *S->qss);
            
            return 1;
        }
        default:
            break;
        }
    }
    return qlua_badconstr(L, "Real");
}

static struct luaL_Reg Qa1(mtLatReal)[] = {
    { "__tostring",   Qa2(q_R,_fmt)      },
    { "__gc",         Qa2(q_R,_gc)       },
    { "__index",      Qa2(q_R,_get)      },
    { "__newindex",   Qa2(q_R,_put)      },
    { "__unm",        Qa2(q_R,_neg)      },
    { "__add",        qlua_add           },
    { "__sub",        qlua_sub           },
    { "__mul",        qlua_mul           },
    { "__div",        qlua_div           },
    { "sum",          Qa2(q_R,_sum)      },
    { "norm2",        Qa2(q_R,_norm2)    },
    { "shift",        Qa2(q_R,_shift)    },
    { "sin",          Qa2(q_R,_sin)      },
    { "cos",          Qa2(q_R,_cos)      },
    { "tan",          Qa2(q_R,_tan)      },
    { "asin",         Qa2(q_R,_asin)     },
    { "acos",         Qa2(q_R,_acos)     },
    { "atan",         Qa2(q_R,_atan)     },
    { "sqrt",         Qa2(q_R,_sqrt)     },
    { "abs",          Qa2(q_R,_abs)      },
    { "exp",          Qa2(q_R,_exp)      },
    { "log",          Qa2(q_R,_log)      },
    { "sign",         Qa2(q_R,_sign)     },
    { "ceil",         Qa2(q_R,_ceil)     },
    { "floor",        Qa2(q_R,_floor)    },
    { "sinh",         Qa2(q_R,_sinh)     },
    { "cosh",         Qa2(q_R,_cosh)     },
    { "tanh",         Qa2(q_R,_tanh)     },
    { "log10",        Qa2(q_R,_log10)    },
    { "expi",         Qa2(q_R,_expi)     },
    { "trunc",        Qa2(q_R,_trunc)    },
    { "round",        Qa2(q_R,_round)    },
    { "set",          Qa2(q_R,_set)      },
    { "float",        Qa2(q_R,_float)    },
    { "double",       Qa2(q_R,_double)   },
    /* "lattice" */
    /* "a-type"  */
    { NULL,           NULL               }
};

Qa1(mLatReal) *
Qa1(qlua_newLatReal)(lua_State *L, int Sidx)
{
    Qa2(QDP_,_Real) *v = Qa2(QDP_,_create_R)();
    Qa1(mLatReal) *hdr;

    if (v == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
        v = Qa2(QDP_,_create_R)();
        if (v == 0)
            luaL_error(L, "not enough memory (QDP_Real)");
    }
    hdr = lua_newuserdata(L, sizeof (Qa1(mLatReal)));
    hdr->ptr = v;
    qlua_createLatticeTable(L, Sidx, Qa1(mtLatReal), Qa1(qLatReal),
                            Qa2(LatReal,Name));
    lua_setmetatable(L, -2);

    return hdr;
}

Qa1(mLatReal) *
Qa1(qlua_checkLatReal)(lua_State *L, int idx, mLattice *S)
{
    void *v = qlua_checkLatticeType(L, idx, Qa1(qLatReal), Qa2(LatReal,Name));
    
    if (S) {
        mLattice *S1 = qlua_ObjLattice(L, idx);
        if (S1->id != S->id)
            luaL_error(L, "%s on a wrong lattice", Qa2(LatReal,Name));
    }

    return (Qa1(mLatReal) *)v;
}

int
Qa2(q_R,_random)(lua_State *L)
{
    mLatRandom *a = qlua_checkLatRandom(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa1(mLatReal) *r = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_R_eq_random_S)(r->ptr, a->ptr, *S->qss);

    return 1;
}

int
Qa2(q_R,_gaussian)(lua_State *L)
{
    mLatRandom *a = qlua_checkLatRandom(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qa1(mLatReal) *r = Qa1(qlua_newLatReal)(L, lua_gettop(L));

    CALL_QDP(L);
    Qa2(QDP_,_R_eq_gaussian_S)(r->ptr, a->ptr, *S->qss);

    return 1;
}

#undef Qab
#undef Qa1
#undef Qa2
#undef Qx
