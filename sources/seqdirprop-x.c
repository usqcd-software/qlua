const char Qs(SeqDirPropName)[] = "qcd.DiracPropagator" Qcolors;
const char Qs(mtnSeqDirProp)[] = "qcd.mtDiracPropagator" Qcolors;

static int
Qs(q_p_fmt)(lua_State *L)
{
    char fmt[72];
    Qs(mSeqDirProp) *b = Qs(qlua_checkSeqDirProp)(L, 1, -1);

    sprintf(fmt, "qcd.DiracPropagator%d(%p)", QC(b), b->ptr);
    lua_pushstring(L, fmt);

    return 1;
}

static int
Qs(q_p_gc)(lua_State *L)
{
    Qs(mSeqDirProp) *b = Qs(qlua_checkSeqDirProp)(L, 1, -1);

    qlua_free(L, b->ptr);
    b->ptr = 0;

    return 0;
}

static int
Qs(q_p_get)(lua_State *L)
{
    switch (qlua_qtype(L, 2)) {
    case qTable: {
        Qs(mSeqDirProp) *V = Qs(qlua_checkSeqDirProp)(L, 1, -1);
        int d = qlua_checkdiracindex(L, 2);
        int c = qlua_checkcolorindex(L, 2, QC(V));
        Qs(mSeqDirFerm) *r = Qs(qlua_newSeqDirFerm)(L, QC(V));
                
        Qx(QLA_D,_D_eq_diracvec_P)(QNC(QC(V)) r->ptr, V->ptr, c, d);
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
Qs(q_p_put)(lua_State *L)
{
    Qs(mSeqDirProp) *V = Qs(qlua_checkSeqDirProp)(L, 1, -1);
    int d = qlua_checkdiracindex(L, 2);
    int c = qlua_checkcolorindex(L, 2, QC(V));
    Qs(mSeqDirFerm) *r = Qs(qlua_checkSeqDirFerm)(L, 3, QC(V));
            
    Qx(QLA_D,_P_eq_diracvec_D)(QNC(QC(V)) V->ptr, r->ptr, c, d);

    return 0;
}

static int
Qs(q_p_norm2_)(lua_State *L)
{
    Qs(mSeqDirProp) *a = Qs(qlua_checkSeqDirProp)(L, 1, -1);
    QLA_D_Real n;

    Qx(QLA_D,_r_eq_norm2_P)(QNC(QC(a)) &n, a->ptr);
    lua_pushnumber(L, n);
    
    return 1;
}

static int
Qs(q_p_conj)(lua_State *L)
{
    Qs(mSeqDirProp) *a = Qs(qlua_checkSeqDirProp)(L, 1, -1);
    Qs(mSeqDirProp) *r = Qs(qlua_newSeqDirProp)(L, QC(a));

    Qx(QLA_D,_P_eq_conj_P)(QNC(QC(a)) r->ptr, a->ptr);
    return 1;
}

static int
Qs(q_p_trans)(lua_State *L)
{
    Qs(mSeqDirProp) *a = Qs(qlua_checkSeqDirProp)(L, 1, -1);
    Qs(mSeqDirProp) *r = Qs(qlua_newSeqDirProp)(L, QC(a));

    Qx(QLA_D,_P_eq_transpose_P)(QNC(QC(a)) r->ptr, a->ptr);
    return 1;
}

static int
Qs(q_p_adjoin)(lua_State *L)
{
    Qs(mSeqDirProp) *a = Qs(qlua_checkSeqDirProp)(L, 1, -1);
    Qs(mSeqDirProp) *r = Qs(qlua_newSeqDirProp)(L, QC(a));

    Qx(QLA_D,_P_eq_Pa)(QNC(QC(a)) r->ptr, a->ptr);
    return 1;
}

static int
Qs(q_p_spintrace)(lua_State *L)
{
    Qs(mSeqDirProp) *a = Qs(qlua_checkSeqDirProp)(L, 1, -1);
    Qs(mSeqColMat) *r = Qs(qlua_newSeqColMat)(L, QC(a));

    Qx(QLA_D,_M_eq_spintrace_P)(QNC(QC(a)) r->ptr, a->ptr);
    return 1;
}

static void
#if QNc == 'N'
Qs(XLA_P_eq_spintranspose_P)(int nc,
                             QLA_DN_DiracPropagator(nc, (*r)),
                             QLA_DN_DiracPropagator(nc, (*a)))
#else
Qs(XLA_P_eq_spintranspose_P)(int nc,
                             Qx(QLA_D,_DiracPropagator) *r,
                             Qx(QLA_D,_DiracPropagator) *a)
#endif
{
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
Qs(q_p_spintranspose)(lua_State *L)
{
    Qs(mSeqDirProp) *a = Qs(qlua_checkSeqDirProp)(L, 1, -1);
    Qs(mSeqDirProp) *r = Qs(qlua_newSeqDirProp)(L, QC(a));

    Qs(XLA_P_eq_spintranspose_P)(QC(a), r->ptr, a->ptr);
    return 1;
}

static void
#if QNc == 'N'
Qs(XLA_C_eq_trace_P)(int nc,
                     QLA_D_Complex *r,
                     QLA_DN_DiracPropagator(nc, (*a)))
#else
Qs(XLA_C_eq_trace_P)(int nc,
                     QLA_D_Complex *r,
                     Qx(QLA_D,_DiracPropagator) *a)
#endif
{
    int is, ic;

    QLA_c_eq_r_plus_ir(*r, 0, 0);
    for (is = 0; is < QDP_Ns; is++) {
        for (ic = 0; ic < nc; ic++)  {
            QLA_c_peq_c(*r, Qx(QLA_D,_elem_P)(*a,ic,is,ic,is));
        }
    }
}

static int
Qs(q_p_trace)(lua_State *L)
{
    Qs(mSeqDirProp) *a = Qs(qlua_checkSeqDirProp)(L, 1, -1);
    QLA_D_Complex *r = qlua_newComplex(L);

    Qs(XLA_C_eq_trace_P)(QC(a), r, a->ptr);
    return 1;
}

static int
Qs(q_p_set)(lua_State *L)
{
    Qs(mSeqDirProp) *r = Qs(qlua_checkSeqDirProp)(L, 1, -1);
    Qs(mSeqDirProp) *a = Qs(qlua_checkSeqDirProp)(L, 2, QC(r));

    Qx(QLA_D,_P_eq_P)(QNC(QC(r)) r->ptr, a->ptr);
    return 0;
}

static int
Qs(q_p_neg)(lua_State *L)
{
    Qs(mSeqDirProp) *a = Qs(qlua_checkSeqDirProp)(L, 1, -1);
    Qs(mSeqDirProp) *r = Qs(qlua_newSeqDirProp)(L, QC(a));

    Qx(QLA_D,_P_eqm_P)(QNC(QC(a)) r->ptr, a->ptr);
    return 1;
}

static int
Qs(q_p_dot_)(lua_State *L)
{
    Qs(mSeqDirProp) *a = Qs(qlua_checkSeqDirProp)(L, 1, -1);
    Qs(mSeqDirProp) *b = Qs(qlua_checkSeqDirProp)(L, 2, QC(a));
    QLA_D_Complex *s = qlua_newComplex(L);

    Qx(QLA_D,_C_eq_P_dot_P)(QNC(QC(a)) s, a->ptr, b->ptr);

    return 1;
}

static int
Qs(q_p_colors)(lua_State *L)
{
#if QNc == 'N'
    Qs(mSeqDirProp) *a = Qs(qlua_checkSeqDirProp)(L, 1, -1);
#else
    Qs(qlua_checkSeqDirProp)(L, 1, -1);
#endif
    lua_pushnumber(L, QC(a));
    return 1;
}

static int
Qs(q_p_copy)(lua_State *L)
{
    Qs(mSeqDirProp) *a = Qs(qlua_checkSeqDirProp)(L, 1, -1);
    Qs(mSeqDirProp) *r = Qs(qlua_newSeqDirProp)(L, QC(a));

    Qx(QLA_D,_P_eq_P)(QNC(QC(a)) r->ptr, a->ptr);
    return 1;
}

static struct luaL_Reg Qs(mtSeqDirProp)[] = {
    { "__tostring",        Qs(q_p_fmt)            },
    { "__gc",              Qs(q_p_gc)             },
    { "__index",           Qs(q_p_get)            },
    { "__newindex",        Qs(q_p_put)            },
    { "__unm",             Qs(q_p_neg)            },
    { "__add",             qlua_add               },
    { "__sub",             qlua_sub               },
    { "__mul",             qlua_mul               },
    { "__div",             qlua_div               },
    { "norm2",             Qs(q_p_norm2_)         },
    { "conj",              Qs(q_p_conj)           },
    { "transpose",         Qs(q_p_trans)          },
    { "adjoin",            Qs(q_p_adjoin)         },
    { "spintrace",         Qs(q_p_spintrace)      },
    { "spintranspose",     Qs(q_p_spintranspose)  },
    { "trace",             Qs(q_p_trace)          },
    { "set",               Qs(q_p_set)            },
    { "colors",            Qs(q_p_colors)         },
    { "copy",              Qs(q_p_copy)           },
    /* "a-type" */
    { NULL,                NULL                   }
};

static int
Qs(q_p_add_p_)(lua_State *L)
{
    Qs(mSeqDirProp) *a = Qs(qlua_checkSeqDirProp)(L, 1, -1);
    Qs(mSeqDirProp) *b = Qs(qlua_checkSeqDirProp)(L, 2, QC(a));
    Qs(mSeqDirProp) *c = Qs(qlua_newSeqDirProp)(L, QC(a));

    Qx(QLA_D,_P_eq_P_plus_P)(QNC(QC(a)) c->ptr, a->ptr, b->ptr);
    return 1;
}

static int
Qs(q_p_sub_p_)(lua_State *L)
{
    Qs(mSeqDirProp) *a = Qs(qlua_checkSeqDirProp)(L, 1, -1);
    Qs(mSeqDirProp) *b = Qs(qlua_checkSeqDirProp)(L, 2, QC(a));
    Qs(mSeqDirProp) *c = Qs(qlua_newSeqDirProp)(L, QC(a));

    Qx(QLA_D,_P_eq_P_minus_P)(QNC(QC(a)) c->ptr, a->ptr, b->ptr);
    return 1;
}

static int
Qs(q_p_mul_r_)(lua_State *L)
{
    Qs(mSeqDirProp) *a = Qs(qlua_checkSeqDirProp)(L, 1, -1);
    QLA_D_Real b = luaL_checknumber(L, 2);
    Qs(mSeqDirProp) *c = Qs(qlua_newSeqDirProp)(L, QC(a));

    Qx(QLA_D,_P_eq_r_times_P)(QNC(QC(a)) c->ptr, &b, a->ptr);
    return 1;
}

static int
Qs(q_r_mul_p_)(lua_State *L)
{
    QLA_D_Real a = luaL_checknumber(L, 1);
    Qs(mSeqDirProp) *b = Qs(qlua_checkSeqDirProp)(L, 2, -1);
    Qs(mSeqDirProp) *c = Qs(qlua_newSeqDirProp)(L, QC(b));

    Qx(QLA_D,_P_eq_r_times_P)(QNC(QC(b)) c->ptr, &a, b->ptr);
    return 1;
}

static int
Qs(q_p_mul_c_)(lua_State *L)
{
    Qs(mSeqDirProp) *a = Qs(qlua_checkSeqDirProp)(L, 1, -1);
    QLA_D_Complex *b = qlua_checkComplex(L, 2);
    Qs(mSeqDirProp) *c = Qs(qlua_newSeqDirProp)(L, QC(a));

    Qx(QLA_D,_P_eq_c_times_P)(QNC(QC(a)) c->ptr, b, a->ptr);
    return 1;
}

static int
Qs(q_c_mul_p_)(lua_State *L)
{
    QLA_D_Complex *a = qlua_checkComplex(L, 1);
    Qs(mSeqDirProp) *b = Qs(qlua_checkSeqDirProp)(L, 2, -1);
    Qs(mSeqDirProp) *c = Qs(qlua_newSeqDirProp)(L, QC(b));

    Qx(QLA_D,_P_eq_c_times_P)(QNC(QC(b)) c->ptr, a, b->ptr);
    return 1;
}

static int
Qs(q_p_mul_p_)(lua_State *L)
{
    Qs(mSeqDirProp) *a = Qs(qlua_checkSeqDirProp)(L, 1, -1);
    Qs(mSeqDirProp) *b = Qs(qlua_checkSeqDirProp)(L, 2, QC(a));
    Qs(mSeqDirProp) *c = Qs(qlua_newSeqDirProp)(L, QC(a));

    Qx(QLA_D,_P_eq_P_times_P)(QNC(QC(a)) c->ptr, a->ptr, b->ptr);
    return 1;
}

static int
Qs(q_p_mul_m_)(lua_State *L)
{
    Qs(mSeqDirProp) *a = Qs(qlua_checkSeqDirProp)(L, 1, -1);
    Qs(mSeqColMat) *b = Qs(qlua_checkSeqColMat)(L, 2, QC(a));
    Qs(mSeqDirProp) *c = Qs(qlua_newSeqDirProp)(L, QC(a));

    Qx(QLA_D,_P_eq_P_times_M)(QNC(QC(a)) c->ptr, a->ptr, b->ptr);
    return 1;
}

static int
Qs(q_m_mul_p_)(lua_State *L)
{
    Qs(mSeqColMat) *a = Qs(qlua_checkSeqColMat)(L, 1, -1);
    Qs(mSeqDirProp) *b = Qs(qlua_checkSeqDirProp)(L, 2, QC(a));
    Qs(mSeqDirProp) *c = Qs(qlua_newSeqDirProp)(L, QC(a));

    Qx(QLA_D,_P_eq_M_times_P)(QNC(QC(a)) c->ptr, a->ptr, b->ptr);
    return 1;
}

static int
Qs(q_p_div_r_)(lua_State *L)
{
    Qs(mSeqDirProp) *a = Qs(qlua_checkSeqDirProp)(L, 1, -1);
    QLA_D_Real b = 1 / luaL_checknumber(L, 2);
    Qs(mSeqDirProp) *c = Qs(qlua_newSeqDirProp)(L, QC(a));

    Qx(QLA_D,_P_eq_r_times_P)(QNC(QC(a)) c->ptr, &b, a->ptr);
    return 1;
}

static int
Qs(q_p_div_c_)(lua_State *L)
{
    Qs(mSeqDirProp) *a = Qs(qlua_checkSeqDirProp)(L, 1, -1);
    QLA_D_Complex *b = qlua_checkComplex(L, 2);
    Qs(mSeqDirProp) *c = Qs(qlua_newSeqDirProp)(L, QC(a));
    double n = 1 / (QLA_real(*b) * QLA_real(*b) + QLA_imag(*b) * QLA_imag(*b));
    QLA_D_Complex s;

    QLA_real(s) = n * QLA_real(*b);
    QLA_imag(s) = -n * QLA_imag(*b);
    Qx(QLA_D,_P_eq_c_times_P)(QNC(QC(a)) c->ptr, &s, a->ptr);
    return 1;
}

Qs(mSeqDirProp) *
Qs(qlua_newSeqDirProp)(lua_State *L, int nc)
{
#if QNc == 'N'
    typedef QLA_DN_DiracPropagator(nc, DataType);
#else
    typedef Qx(QLA_D,_DiracPropagator) DataType;
#endif
    Qs(mSeqDirProp) *hdr = lua_newuserdata(L, sizeof (Qs(mSeqDirFerm)));
    hdr->ptr = qlua_malloc(L, sizeof (DataType));
#if QNc == 'N'
    hdr->nc = nc;
#endif
    luaL_getmetatable(L, Qs(mtnSeqDirProp));
    lua_setmetatable(L, -2);

    return hdr;
}

Qs(mSeqDirProp) *
Qs(qlua_newZeroSeqDirProp)(lua_State *L, int nc)
{
        Qs(mSeqDirProp) *v = Qs(qlua_newSeqDirProp)(L, nc);
        Qx(QLA_D,_P_eq_zero)(QNC(nc) v->ptr);
        return v;
}

Qs(mSeqDirProp) *
Qs(qlua_checkSeqDirProp)(lua_State *L, int idx, int nc)
{
    void *v = qlua_checkLatticeType(L, idx, Qs(qSeqDirProp),
                                    Qs(SeqDirPropName));
    Qs(mSeqDirProp) *z = (Qs(mSeqDirProp) *)v;
#if QNc == 'N'
    if (nc != -1) {
        if (z->nc != nc)
            luaL_error(L, "Wrong number of colors");
    }
#endif

    return z;
}

static int
Qs(q_seqdirprop_)(lua_State *L, int nc)
{
    switch (lua_gettop(L)) {
    case 1: {
        Qs(qlua_newZeroSeqDirProp)(L, nc);
        return 1;
    }
    case 2: {
        Qs(mSeqDirProp) *v = Qs(qlua_newZeroSeqDirProp)(L, nc);
        
        switch (qlua_qtype(L, 2)) {
        case Qs(qSeqColMat): {
            Qs(mSeqColMat) *w = Qs(qlua_checkSeqColMat)(L, 2, nc);
            QLA_D_Complex c;
            int ic, jc, ks;

            for (ic = 0; ic < nc; ic++) {
                for (jc = 0; jc < nc; jc++) {
                    Qx(QLA_D,_C_eq_elem_M)(QNC(nc) &c, w->ptr, ic, jc);
                    for (ks = 0; ks < QDP_Ns; ks++)
                        Qx(QLA_D,_P_eq_elem_C)(QNC(nc) v->ptr, &c,
                                               ic, ks, jc, ks);
                }
            }
            lua_pop(L, 1);

            return 1;
        }
        case qReal: {
          double r = luaL_checknumber(L, 2);
          QLA_D_Complex c;
          int ic, is;
          
          QLA_real(c) = r;
          QLA_imag(c) = 0.0;
          for (ic = 0; ic < nc; ic++) {
            for (is = 0; is < QDP_Ns; is++) {
              Qx(QLA_D,_P_eq_elem_C)(QNC(nc) v->ptr, &c, ic, is, ic, is);
              }
          }
          return 1;
        }
        case qComplex: {
          QLA_D_Complex *c = qlua_checkComplex(L, 2);
          int ic, is;
          
          for (ic = 0; ic < nc; ic++) {
            for (is = 0; is < QDP_Ns; is++) {
              Qx(QLA_D,_P_eq_elem_C)(QNC(nc) v->ptr, c, ic, is, ic, is);
            }
          }
          return 1;
        }
        default:
            break;
        }
    }
    case 3: {
        Qs(mSeqDirFerm) *z = Qs(qlua_checkSeqDirFerm)(L, 2, nc);
        int d = qlua_checkdiracindex(L, 3);
        int c = qlua_checkcolorindex(L, 3, nc);
        Qs(mSeqDirProp) *v = Qs(qlua_newZeroSeqDirProp)(L, nc);

        Qx(QLA_D,_P_eq_diracvec_D)(QNC(nc) v->ptr, z->ptr, c, d);
        return 1;
    }
    }
    return qlua_badconstr(L, "DiracPropagator" Qcolors);
}

static const QLUA_Op2 Qs(ops)[] = {
    { qlua_add_table, Qs(qSeqDirProp),  Qs(qSeqDirProp),  Qs(q_p_add_p_) },
    { qlua_sub_table, Qs(qSeqDirProp),  Qs(qSeqDirProp),  Qs(q_p_sub_p_) },
    { qlua_mul_table, qReal,            Qs(qSeqDirProp),  Qs(q_r_mul_p_) },
    { qlua_mul_table, Qs(qSeqDirProp),  qReal,            Qs(q_p_mul_r_) },
    { qlua_mul_table, qComplex,         Qs(qSeqDirProp),  Qs(q_c_mul_p_) },
    { qlua_mul_table, Qs(qSeqDirProp),  qComplex,         Qs(q_p_mul_c_) },
    { qlua_mul_table, Qs(qSeqDirProp),  Qs(qSeqDirProp),  Qs(q_p_mul_p_) },
    { qlua_mul_table, Qs(qSeqDirProp),  Qs(qSeqColMat),   Qs(q_p_mul_m_) },
    { qlua_mul_table, Qs(qSeqColMat),   Qs(qSeqDirProp),  Qs(q_m_mul_p_) },
    { qlua_div_table, Qs(qSeqDirProp),  qReal,            Qs(q_p_div_r_) },
    { qlua_div_table, Qs(qSeqDirProp),  qComplex,         Qs(q_p_div_c_) },
    { NULL,           qNoType,          qNoType,          NULL           }
};

#undef QNc
#undef Qcolors
#undef Qs
#undef Qx
#undef QC
#undef QNC
