static const char Qs(SeqDirFermName)[] = "qcd.DiracFermion" Qcolors;
static const char Qs(mtnSeqDirFerm)[] = "qcd.mtDiracFermion" Qcolors;

static int
Qs(q_d_fmt)(lua_State *L)
{
    char fmt[72];
    Qs(mSeqDirFerm) *b = Qs(qlua_checkSeqDirFerm)(L, 1, -1);

    sprintf(fmt, "qcd.DiracFermion%d(%p)", QC(b), b->ptr);
    lua_pushstring(L, fmt);

    return 1;
}

static int
Qs(q_d_gc)(lua_State *L)
{
    Qs(mSeqDirFerm) *b = Qs(qlua_checkSeqDirFerm)(L, 1, -1);

    qlua_free(L, b->ptr);
    b->ptr = 0;

    return 0;
}

static int
Qs(q_d_get)(lua_State *L)
{
    switch (qlua_qtype(L, 2)) {
    case qTable: {
        Qs(mSeqDirFerm) *V = Qs(qlua_checkSeqDirFerm)(L, 1, -1);
        int d = qlua_checkdiracindex(L, 2);
        int c = qlua_colorindex(L, 2, QC(V));
        if (c == -1) {
            Qs(mSeqColVec) *r = Qs(qlua_newSeqColVec)(L, QC(V));

            Qx(QLA_D,_V_eq_colorvec_D)(QNC(QC(V)) r->ptr, V->ptr, d);
        } else {
            QLA_D_Complex *r = qlua_newComplex(L);

            Qx(QLA_D,_C_eq_elem_D)(QNC(QC(V)) r, V->ptr, c, d);
        }
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
Qs(q_d_put)(lua_State *L)
{
    Qs(mSeqDirFerm) *V = Qs(qlua_checkSeqDirFerm)(L, 1, -1);
    int d = qlua_checkdiracindex(L, 2);
    int c = qlua_colorindex(L, 2, QC(V));

    if (c == -1) {
        Qs(mSeqColVec) *z = Qs(qlua_checkSeqColVec)(L, 3, QC(V));
        
        Qx(QLA_D,_D_eq_colorvec_V)(QNC(QC(V)) V->ptr, z->ptr, d);
    } else {
        QLA_D_Complex *z = qlua_checkComplex(L, 3);
        
        Qx(QLA_D,_D_eq_elem_C)(QNC(QC(V)) V->ptr, z, c, d);
    }
    return 0;
}

static int
Qs(q_d_neg)(lua_State *L)
{
    Qs(mSeqDirFerm) *a = Qs(qlua_checkSeqDirFerm)(L, 1, -1);
    Qs(mSeqDirFerm) *r = Qs(qlua_newSeqDirFerm)(L, QC(a));
    QLA_D_Real m1 = -1;

    Qx(QLA_D,_D_eq_r_times_D)(QNC(QC(a)) r->ptr, &m1, a->ptr);
    return 1;
}

static int
Qs(q_d_norm2_)(lua_State *L)
{
    Qs(mSeqDirFerm) *a = Qs(qlua_checkSeqDirFerm)(L, 1, -1);
    QLA_D_Real n;

    Qx(QLA_D,_r_eq_norm2_D)(QNC(QC(a)) &n, a->ptr);
    lua_pushnumber(L, n);
    return 1;
}

static int
Qs(q_d_conj)(lua_State *L)
{
    Qs(mSeqDirFerm) *a = Qs(qlua_checkSeqDirFerm)(L, 1, -1);
    Qs(mSeqDirFerm) *r = Qs(qlua_newSeqDirFerm)(L, QC(a));

    Qx(QLA_D,_D_eq_conj_D)(QNC(QC(a)) r->ptr, a->ptr);
    return 1;
}

static int
Qs(q_d_gamma)(lua_State *L)
{
    Qs(mSeqDirFerm) *f = Qs(qlua_checkSeqDirFerm)(L, 1, -1);
    int mu = qlua_gammaindex(L, 2);
    int n = qlua_gammabinary(L, 2);
    
    if (((n == -1) && (mu == -1)) || ((n != -1) && (mu != -1)))
        return qlua_badindex(L, "DiracFermion" Qcolors);
    if (n == -1) {
        Qs(mSeqDirFerm) *r = Qs(qlua_newSeqDirFerm)(L, QC(f));

        if (mu < 5) {
            Qx(QLA_D,_D_eq_gamma_times_D)(QNC(QC(f)) r->ptr, f->ptr, 1 << mu);
        } else {
            Qx(QLA_D,_D_eq_gamma_times_D)(QNC(QC(f)) r->ptr, f->ptr, 15);
        }
        return 1;
    }
    if (mu == -1) {
        Qs(mSeqDirFerm) *r = Qs(qlua_newSeqDirFerm)(L, QC(f));
        
        Qx(QLA_D,_D_eq_gamma_times_D)(QNC(QC(f)) r->ptr, f->ptr, n);
        return 1;
    }

    return luaL_error(L, "this could never happen");
}

static int
Qs(q_d_set)(lua_State *L)
{
    Qs(mSeqDirFerm) *r = Qs(qlua_checkSeqDirFerm)(L, 1, -1);
    Qs(mSeqDirFerm) *a = Qs(qlua_checkSeqDirFerm)(L, 2, QC(r));

    Qx(QLA_D,_D_eq_D)(QNC(QC(r)) r->ptr, a->ptr);
    lua_pop(L, 2);

    return 0;
}

static int
Qs(q_d_dot_)(lua_State *L)
{
    Qs(mSeqDirFerm) *a = Qs(qlua_checkSeqDirFerm)(L, 1, -1);
    Qs(mSeqDirFerm) *b = Qs(qlua_checkSeqDirFerm)(L, 2, QC(a));
    QLA_D_Complex *s = qlua_newComplex(L);

    Qx(QLA_D,_C_eq_D_dot_D)(QNC(QC(a)) s, a->ptr, b->ptr);

    return 1;
}

static int
Qs(q_d_add_d_)(lua_State *L)
{
    Qs(mSeqDirFerm) *a = Qs(qlua_checkSeqDirFerm)(L, 1, -1);
    Qs(mSeqDirFerm) *b = Qs(qlua_checkSeqDirFerm)(L, 2, QC(a));
    Qs(mSeqDirFerm) *c = Qs(qlua_newSeqDirFerm)(L, QC(a));

    Qx(QLA_D,_D_eq_D_plus_D)(QNC(QC(a)) c->ptr, a->ptr, b->ptr);
    return 1;
}

static int
Qs(q_d_sub_d_)(lua_State *L)
{
    Qs(mSeqDirFerm) *a = Qs(qlua_checkSeqDirFerm)(L, 1, -1);
    Qs(mSeqDirFerm) *b = Qs(qlua_checkSeqDirFerm)(L, 2, QC(a));
    Qs(mSeqDirFerm) *c = Qs(qlua_newSeqDirFerm)(L, QC(a));

    Qx(QLA_D,_D_eq_D_minus_D)(QNC(QC(a)) c->ptr, a->ptr, b->ptr);
    return 1;
}

static int
Qs(q_d_mul_r_)(lua_State *L)
{
    Qs(mSeqDirFerm) *a = Qs(qlua_checkSeqDirFerm)(L, 1, -1);
    QLA_D_Real b = luaL_checknumber(L, 2);
    Qs(mSeqDirFerm) *c = Qs(qlua_newSeqDirFerm)(L, QC(a));

    Qx(QLA_D,_D_eq_r_times_D)(QNC(QC(a)) c->ptr, &b, a->ptr);
    return 1;
}

static int
Qs(q_r_mul_d_)(lua_State *L)
{
    QLA_D_Real a = luaL_checknumber(L, 1);
    Qs(mSeqDirFerm) *b = Qs(qlua_checkSeqDirFerm)(L, 2, -1);
    Qs(mSeqDirFerm) *c = Qs(qlua_newSeqDirFerm)(L, QC(b));

    Qx(QLA_D,_D_eq_r_times_D)(QNC(QC(b)) c->ptr, &a, b->ptr);
    return 1;
}

static int
Qs(q_d_mul_c_)(lua_State *L)
{
    Qs(mSeqDirFerm) *a = Qs(qlua_checkSeqDirFerm)(L, 1, -1);
    QLA_D_Complex *b = qlua_checkComplex(L, 2);
    Qs(mSeqDirFerm) *c = Qs(qlua_newSeqDirFerm)(L, QC(a));

    Qx(QLA_D,_D_eq_c_times_D)(QNC(QC(a)) c->ptr, b, a->ptr);
    return 1;
}

static int
Qs(q_c_mul_d_)(lua_State *L)
{
    QLA_D_Complex *a = qlua_checkComplex(L, 1);
    Qs(mSeqDirFerm) *b = Qs(qlua_checkSeqDirFerm)(L, 2, -1);
    Qs(mSeqDirFerm) *c = Qs(qlua_newSeqDirFerm)(L, QC(b));

    Qx(QLA_D,_D_eq_c_times_D)(QNC(QC(b)) c->ptr, a, b->ptr);
    return 1;
}

static int
Qs(q_m_mul_d_)(lua_State *L)
{
    Qs(mSeqColMat) *a = Qs(qlua_checkSeqColMat)(L, 1, -1);
    Qs(mSeqDirFerm) *b = Qs(qlua_checkSeqDirFerm)(L, 2, QC(a));
    Qs(mSeqDirFerm) *c = Qs(qlua_newSeqDirFerm)(L, QC(a));

    Qx(QLA_D,_D_eq_M_times_D)(QNC(QC(a)) c->ptr, a->ptr, b->ptr);
    return 1;
}

static int
Qs(q_d_div_r_)(lua_State *L)
{
    Qs(mSeqDirFerm) *a = Qs(qlua_checkSeqDirFerm)(L, 1, -1);
    QLA_D_Real b = 1 / luaL_checknumber(L, 2);
    Qs(mSeqDirFerm) *c = Qs(qlua_newSeqDirFerm)(L, QC(a));

    Qx(QLA_D,_D_eq_r_times_D)(QNC(QC(a)) c->ptr, &b, a->ptr);
    return 1;
}

static int
Qs(q_d_div_c_)(lua_State *L)
{
    Qs(mSeqDirFerm) *a = Qs(qlua_checkSeqDirFerm)(L, 1, -1);
    QLA_D_Complex *b = qlua_checkComplex(L, 2);
    Qs(mSeqDirFerm) *c = Qs(qlua_newSeqDirFerm)(L, QC(a));
    double n = 1 / (QLA_real(*b) * QLA_real(*b) + QLA_imag(*b) * QLA_imag(*b));
    QLA_D_Complex s;

    QLA_real(s) = n * QLA_real(*b);
    QLA_imag(s) = -n * QLA_imag(*b);
    Qx(QLA_D,_D_eq_c_times_D)(QNC(QC(a)) c->ptr, &s, a->ptr);
    return 1;
}

static int
Qs(q_d_colors)(lua_State *L)
{
#if QNc == 'N'
    Qs(mSeqDirFerm) *a = Qs(qlua_checkSeqDirFerm)(L, 1, -1);
#else
    Qs(qlua_checkSeqDirFerm)(L, 1, -1);
#endif
    lua_pushnumber(L, QC(a));

    return 1;
}

static int
Qs(q_d_copy)(lua_State *L)
{
    Qs(mSeqDirFerm) *a = Qs(qlua_checkSeqDirFerm)(L, 1, -1);
    Qs(mSeqDirFerm) *r = Qs(qlua_newSeqDirFerm)(L, QC(a));

    Qx(QLA_D,_D_eq_D)(QNC(QC(a)) r->ptr, a->ptr);
    return 1;
}

static struct luaL_Reg Qs(mtSeqDirFerm)[] = {
    { "__tostring",        Qs(q_d_fmt)     },
    { "__gc",              Qs(q_d_gc)      },
    { "__index",           Qs(q_d_get)     },
    { "__newindex",        Qs(q_d_put)     },
    { "__unm",             Qs(q_d_neg)     },
    { "__add",             qlua_add        },
    { "__sub",             qlua_sub        },
    { "__mul",             qlua_mul        },
    { "__div",             qlua_div        },
    { "norm2",             Qs(q_d_norm2_)  },
    { "conj",              Qs(q_d_conj)    },
    { "gamma",             Qs(q_d_gamma)   },
    { "set",               Qs(q_d_set)     },
    { "colors",            Qs(q_d_colors)  },
    { "copy",              Qs(q_d_copy)    },
    /* "a-type" */
    { NULL,                NULL }
};

Qs(mSeqDirFerm) *
Qs(qlua_newSeqDirFerm)(lua_State *L, int nc)
{
#if QNc == 'N'
    typedef QLA_DN_DiracFermion(nc, DataType);
#else
    typedef Qx(QLA_D,_DiracFermion) DataType;
#endif
    Qs(mSeqDirFerm) *hdr = lua_newuserdata(L, sizeof (Qs(mSeqDirFerm)));
    hdr->ptr = qlua_malloc(L, sizeof (DataType));
#if QNc == 'N'
    hdr->nc = nc;
#endif
    luaL_getmetatable(L, Qs(mtnSeqDirFerm));
    lua_setmetatable(L, -2);

    return hdr;
}

Qs(mSeqDirFerm) *
Qs(qlua_checkSeqDirFerm)(lua_State *L, int idx, int nc)
{
    void *v = qlua_checkLatticeType(L, idx, Qs(qSeqDirFerm),
                                    Qs(SeqDirFermName));
    Qs(mSeqDirFerm) *z = (Qs(mSeqDirFerm) *)v;
#if QNc == 'N'
    if (nc != -1) {
        if (z->nc != nc)
            luaL_error(L, "Wrong number of colors");
    }
#endif

    return z;
}

static int
Qs(q_seqdirferm_)(lua_State *L, int nc)
{
    switch (lua_gettop(L)) {
    case 0: {
        Qs(mSeqDirFerm) *v = Qs(qlua_newSeqDirFerm)(L, nc);

        Qx(QLA_D,_D_eq_zero)(QNC(nc) v->ptr);
        return 1;
    }
    case 2: {
        switch (qlua_qtype(L, 1)) {
        case qComplex: {
            QLA_D_Complex *z = qlua_checkComplex(L, 1);
            int c = qlua_checkcolorindex(L, 2, nc);
            int d = qlua_checkdiracindex(L, 2);
            Qs(mSeqDirFerm) *v = Qs(qlua_newSeqDirFerm)(L, nc);

            Qx(QLA_D,_D_eq_zero)(QNC(nc) v->ptr);
            Qx(QLA_D,_D_eq_elem_C)(QNC(nc) v->ptr, z, c, d);

            return 1;
        }
        case Qs(qSeqColVec): {
            Qs(mSeqColVec) *w = Qs(qlua_checkSeqColVec)(L, 1, nc);
            int d = qlua_checkdiracindex(L, 2);
            Qs(mSeqDirFerm) *v = Qs(qlua_newSeqDirFerm)(L, nc);

            Qx(QLA_D,_D_eq_zero)(QNC(nc) v->ptr);
            Qx(QLA_D,_D_eq_colorvec_V)(QNC(nc) v->ptr, w->ptr, d);
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
    { qlua_add_table, Qs(qSeqDirFerm),  Qs(qSeqDirFerm),  Qs(q_d_add_d_) },
    { qlua_sub_table, Qs(qSeqDirFerm),  Qs(qSeqDirFerm),  Qs(q_d_sub_d_) },
    { qlua_mul_table, qReal,            Qs(qSeqDirFerm),  Qs(q_r_mul_d_) },
    { qlua_mul_table, Qs(qSeqDirFerm),  qReal,            Qs(q_d_mul_r_) },
    { qlua_mul_table, qComplex,         Qs(qSeqDirFerm),  Qs(q_c_mul_d_) },
    { qlua_mul_table, Qs(qSeqDirFerm),  qComplex,         Qs(q_d_mul_c_) },
    { qlua_mul_table, Qs(qSeqColMat),   Qs(qSeqDirFerm),  Qs(q_m_mul_d_) },
    { qlua_div_table, Qs(qSeqDirFerm),  qReal,            Qs(q_d_div_r_) },
    { qlua_div_table, Qs(qSeqDirFerm),  qComplex,         Qs(q_d_div_c_) },
    { NULL,           qNoType,          qNoType,          NULL           }
};

#undef QNc
#undef Qcolors
#undef Qs
#undef Qx
#undef QC
#undef QNC
