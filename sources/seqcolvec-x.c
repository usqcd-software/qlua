static const char Qs(SeqColVecName)[] = "qcd.ColorVector" Qcolors;
static const char Qs(mtnSeqColVec)[] = "qcd.mtColorVector" Qcolors;

static int
Qs(q_v_fmt)(lua_State *L)
{
    char fmt[72];
    Qs(mSeqColVec) *b = Qs(qlua_checkSeqColVec)(L, 1, -1);

    sprintf(fmt, "qcd.ColorVector%d(%p)", QC(b), b->ptr);
    lua_pushstring(L, fmt);

    return 1;
}

static int
Qs(q_v_gc)(lua_State *L)
{
    Qs(mSeqColVec) *b = Qs(qlua_checkSeqColVec)(L, 1, -1);

    qlua_free(L, b->ptr);
    b->ptr = 0;

    return 0;
}

static int
Qs(q_v_get)(lua_State *L)
{
    switch (qlua_qtype(L, 2)) {
    case qTable: {
        Qs(mSeqColVec) *V = Qs(qlua_checkSeqColVec)(L, 1, -1);
        int c = qlua_checkcolorindex(L, 2, QC(V));
        QLA_D_Complex *r = qlua_newComplex(L);

        Qx(QLA_D,_C_eq_elem_V)(QNC(QC(V)) r, V->ptr, c);
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
Qs(q_v_put)(lua_State *L)
{
    Qs(mSeqColVec) *V = Qs(qlua_checkSeqColVec)(L, 1, -1);
    int c = qlua_checkcolorindex(L, 2, QC(V));
    QLA_D_Complex *z = qlua_checkComplex(L, 3);
        
    Qx(QLA_D,_V_eq_elem_C)(QNC(QC(V)) V->ptr, z, c);

    return 0;
}

static int
Qs(q_v_dot_)(lua_State *L)
{
    Qs(mSeqColVec) *a = Qs(qlua_checkSeqColVec)(L, 1, -1);
    Qs(mSeqColVec) *b = Qs(qlua_checkSeqColVec)(L, 2, QC(a));
    QLA_D_Complex *s = qlua_newComplex(L);

#if QNc == 'N'
    Qx(QLA_D,_C_eq_V_dot_V)(QC(a), s, a->ptr, b->ptr);
#else
    Qx(QLA_D,_C_eq_V_dot_V)(s, a->ptr, b->ptr);
#endif

    return 1;
}

static int
Qs(q_v_add_v_)(lua_State *L)
{
    Qs(mSeqColVec) *a = Qs(qlua_checkSeqColVec)(L, 1, -1);
    Qs(mSeqColVec) *b = Qs(qlua_checkSeqColVec)(L, 2, QC(a));
    Qs(mSeqColVec) *c = Qs(qlua_newSeqColVec)(L, QC(a));

    Qx(QLA_D,_V_eq_V_plus_V)(QNC(QC(a)) c->ptr, a->ptr, b->ptr);

    return 1;
}

static int
Qs(q_v_sub_v_)(lua_State *L)
{
    Qs(mSeqColVec) *a = Qs(qlua_checkSeqColVec)(L, 1, -1);
    Qs(mSeqColVec) *b = Qs(qlua_checkSeqColVec)(L, 2, QC(a));
    Qs(mSeqColVec) *c = Qs(qlua_newSeqColVec)(L, QC(a));

    Qx(QLA_D,_V_eq_V_minus_V)(QNC(QC(a)) c->ptr, a->ptr, b->ptr);

    return 1;
}

static int
Qs(q_v_mul_r_)(lua_State *L)
{
    Qs(mSeqColVec) *a = Qs(qlua_checkSeqColVec)(L, 1, -1);
    QLA_Real b = luaL_checknumber(L, 2);
    Qs(mSeqColVec) *c = Qs(qlua_newSeqColVec)(L, QC(a));

    Qx(QLA_D,_V_eq_r_times_V)(QNC(QC(a)) c->ptr, &b, a->ptr);

    return 1;
}

static int
Qs(q_r_mul_v_)(lua_State *L)
{
    QLA_Real a = luaL_checknumber(L, 1);
    Qs(mSeqColVec) *b = Qs(qlua_checkSeqColVec)(L, 2, -1);
    Qs(mSeqColVec) *c = Qs(qlua_newSeqColVec)(L, QC(b));

    Qx(QLA_D,_V_eq_r_times_V)(QNC(QC(b)) c->ptr, &a, b->ptr);

    return 1;
}

static int
Qs(q_v_mul_c_)(lua_State *L)
{
    Qs(mSeqColVec) *a = Qs(qlua_checkSeqColVec)(L, 1, -1);
    QLA_Complex *b = qlua_checkComplex(L, 2);
    Qs(mSeqColVec) *c = Qs(qlua_newSeqColVec)(L, QC(a));

    Qx(QLA_D,_V_eq_c_times_V)(QNC(QC(a)) c->ptr, b, a->ptr);

    return 1;
}

static int
Qs(q_c_mul_v_)(lua_State *L)
{
    QLA_Complex *a = qlua_checkComplex(L, 1);
    Qs(mSeqColVec) *b = Qs(qlua_checkSeqColVec)(L, 2, -1);
    Qs(mSeqColVec) *c = Qs(qlua_newSeqColVec)(L, QC(b));

    Qx(QLA_D,_V_eq_c_times_V)(QNC(QC(b)) c->ptr, a, b->ptr);

    return 1;
}

static int
Qs(q_v_norm2_)(lua_State *L)
{
    Qs(mSeqColVec) *a = Qs(qlua_checkSeqColVec)(L, 1, -1);
    QLA_D_Real n;

    Qx(QLA_D,_r_eq_norm2_V)(QNC(QC(a)) &n, a->ptr);
    lua_pushnumber(L, n);
    
    return 1;
}

static int
Qs(q_v_conj)(lua_State *L)
{
    Qs(mSeqColVec) *a = Qs(qlua_checkSeqColVec)(L, 1, -1);
    Qs(mSeqColVec) *r = Qs(qlua_newSeqColVec)(L, QC(a));

    Qx(QLA_D,_V_eq_conj_V)(QNC(QC(a)) r->ptr, a->ptr);

    return 1;
}

static int
Qs(q_v_set)(lua_State *L)
{
    Qs(mSeqColVec) *r = Qs(qlua_checkSeqColVec)(L, 1, -1);
    Qs(mSeqColVec) *a = Qs(qlua_checkSeqColVec)(L, 2, QC(r));

    Qx(QLA_D,_V_eq_V)(QNC(QC(r)) r->ptr, a->ptr);
    lua_pop(L, 2);

    return 0;
}

static int
Qs(q_v_neg)(lua_State *L)
{
    Qs(mSeqColVec) *a = Qs(qlua_checkSeqColVec)(L, 1, -1);
    Qs(mSeqColVec) *r = Qs(qlua_newSeqColVec)(L, QC(a));
    QLA_D_Real m1 = -1;

    Qx(QLA_D,_V_eq_r_times_V)(QNC(QC(a)) r->ptr, &m1, a->ptr);

    return 1;
}

static int
Qs(q_v_div_r_)(lua_State *L)
{
    Qs(mSeqColVec) *a = Qs(qlua_checkSeqColVec)(L, 1, -1);
    QLA_D_Real b = 1 / luaL_checknumber(L, 2);
    Qs(mSeqColVec) *c = Qs(qlua_newSeqColVec)(L, QC(a));

    Qx(QLA_D,_V_eq_r_times_V)(QNC(QC(a)) c->ptr, &b, a->ptr);

    return 1;
}

static int
Qs(q_v_div_c_)(lua_State *L)
{
    Qs(mSeqColVec) *a = Qs(qlua_checkSeqColVec)(L, 1, -1);
    QLA_D_Complex *b = qlua_checkComplex(L, 2);
    Qs(mSeqColVec) *c = Qs(qlua_newSeqColVec)(L, QC(a));
    double n = 1 / (QLA_real(*b) * QLA_real(*b) + QLA_imag(*b) * QLA_imag(*b));
    QLA_D_Complex s;

    QLA_real(s) = n * QLA_real(*b);
    QLA_imag(s) = -n * QLA_imag(*b);
    Qx(QLA_D,_V_eq_c_times_V)(QNC(QC(a)) c->ptr, &s, a->ptr);

    return 1;
}

static int
Qs(q_v_colors)(lua_State *L)
{
#if QNc == 'N'
    Qs(mSeqColVec) *a = Qs(qlua_checkSeqColVec)(L, 1, -1);
#else
    Qs(qlua_checkSeqColVec)(L, 1, -1);
#endif
    lua_pushnumber(L, QC(a));

    return 1;
}

static int
Qs(q_v_copy)(lua_State *L)
{
    Qs(mSeqColVec) *a = Qs(qlua_checkSeqColVec)(L, 1, -1);
    Qs(mSeqColVec) *r = Qs(qlua_newSeqColVec)(L, QC(a));

    Qx(QLA_D,_V_eq_V)(QNC(QC(a)) r->ptr, a->ptr);

    return 1;
}
static struct luaL_Reg Qs(mtSeqColVec)[] = {
    { "__tostring",        Qs(q_v_fmt)    },
    { "__gc",              Qs(q_v_gc)     },
    { "__index",           Qs(q_v_get)    },
    { "__newindex",        Qs(q_v_put)    },
    { "__unm",             Qs(q_v_neg)    },
    { "__add",             qlua_add       },
    { "__sub",             qlua_sub       },
    { "__mul",             qlua_mul       },
    { "__div",             qlua_div       },
    { "norm2",             Qs(q_v_norm2_) },
    { "conj",              Qs(q_v_conj)   },
    { "set",               Qs(q_v_set)    },
    { "colors",            Qs(q_v_colors) },
    { "copy",              Qs(q_v_copy)   },
    /* "a-type" */
    { NULL,                NULL           }
};

Qs(mSeqColVec) *
Qs(qlua_newSeqColVec)(lua_State *L, int nc)
{
#if QNc == 'N'
    typedef QLA_DN_ColorVector(nc, DataType);
#else
    typedef Qx(QLA_D,_ColorVector) DataType;
#endif
    Qs(mSeqColVec) *hdr = lua_newuserdata(L, sizeof (Qs(mSeqColVec)));
    hdr->ptr = qlua_malloc(L, sizeof (DataType));
#if QNc == 'N'
    hdr->nc = nc;
#endif
    luaL_getmetatable(L, Qs(mtnSeqColVec));
    lua_setmetatable(L, -2);

    return hdr;
}

Qs(mSeqColVec) *
Qs(qlua_newZeroSeqColVec)(lua_State *L, int nc)
{
        Qs(mSeqColVec) *v = Qs(qlua_newSeqColVec)(L, nc);
        Qx(QLA_D,_V_eq_zero)(QNC(nc) v->ptr);
        return v;
}
Qs(mSeqColVec) *
Qs(qlua_checkSeqColVec)(lua_State *L, int idx, int nc)
{
    void *v = qlua_checkLatticeType(L, idx, Qs(qSeqColVec), Qs(SeqColVecName));
    Qs(mSeqColVec) *z = (Qs(mSeqColVec) *)v;

#if QNc == 'N'
    if (nc != -1) {
        if (z->nc != nc)
            luaL_error(L, "Wrong number of colors");
    }
#endif
    return z;
}

static int
Qs(q_seqcolvec_)(lua_State *L, int nc)
{
    switch (lua_gettop(L)) {
    case 1: {
        Qs(qlua_newZeroSeqColVec)(L, nc);
        return 1;
    }
    case 2: {
        Qs(mSeqColVec) *a = Qs(qlua_checkSeqColVec)(L, 2, nc);
        Qs(mSeqColVec) *r = Qs(qlua_newSeqColVec)(L, nc);
        
        Qx(QLA_D,_V_eq_V)(QNC(nc) r->ptr, a->ptr);
        
        return 1;
    }
    case 3: {
        QLA_D_Complex *c = qlua_checkComplex(L, 2);
        int a = qlua_checkcolorindex(L, 3, nc);
        Qs(mSeqColVec) *r = Qs(qlua_newSeqColVec)(L, nc);

        Qx(QLA_D,_V_eq_elem_C)(QNC(nc) r->ptr, c, a);

        return 1;
    }
    }
    return qlua_badconstr(L, "ColorVector" Qcolors);
}

static const QLUA_Op2 Qs(ops)[] = {
    { qlua_add_table, Qs(qSeqColVec),  Qs(qSeqColVec),  Qs(q_v_add_v_) },
    { qlua_sub_table, Qs(qSeqColVec),  Qs(qSeqColVec),  Qs(q_v_sub_v_) },
    { qlua_mul_table, qReal,           Qs(qSeqColVec),  Qs(q_r_mul_v_) },
    { qlua_mul_table, Qs(qSeqColVec),  qReal,           Qs(q_v_mul_r_) },
    { qlua_mul_table, qComplex,        Qs(qSeqColVec),  Qs(q_c_mul_v_) },
    { qlua_mul_table, Qs(qSeqColVec),  qComplex,        Qs(q_v_mul_c_) },
    { qlua_div_table, Qs(qSeqColVec),  qReal,           Qs(q_v_div_r_) },
    { qlua_div_table, Qs(qSeqColVec),  qComplex,        Qs(q_v_div_c_) },
    { NULL,           qNoType,         qNoType,         NULL        }
    
};
#undef QNc
#undef Qcolors
#undef Qs
#undef Qx
#undef QC
#undef QNC
