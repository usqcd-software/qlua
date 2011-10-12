
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

static void
#if QNc == 'N'
Qs(unpack_vec)(int nc, double *std, QLA_DN_ColorVector(nc, (*src)))
#else
Qs(unpack_vec)(int nc, double *std, Qx(QLA_D,_ColorVector) *src)
#endif
{
    int c;

    for (c = 0; c < nc; c++) {
        QLA_D_Complex *z = &Qx(QLA_D,_elem_V)(*src, c);
        std[2*c] = QLA_real(*z);
        std[2*c + 1] = QLA_imag(*z);
    }
}

static void
#if QNc == 'N'
Qs(pack_vec)(int nc, QLA_DN_ColorVector(nc, (*dst)), const double *std)
#else
Qs(pack_vec)(int nc, Qx(QLA_D,_ColorVector) *dst, const double *std)
#endif
{
    QLA_D_Complex z;
    int c;

    for (c = 0; c < nc; c++) {
        QLA_c_eq_r_plus_ir(z, std[2 * c], std[2 * c + 1]);
        QLA_c_eq_c(Qx(QLA_D,_elem_V)(*dst, c), z);
    }
}

static int
Qs(q_V_get)(lua_State *L)
{
    switch (qlua_qtype(L, 2)) {
    case qTable: {
        Qs(mLatColVec) *V = Qs(qlua_checkLatColVec)(L, 1, NULL, -1);
        mLattice *S = qlua_ObjLattice(L, 1);
        int *idx = qlua_latcoord(L, 2, S);
        int Sidx = lua_gettop(L);

        if (idx == NULL) {
            /* V[{c=val}] is the only acceptable form here */
            int c = qlua_checkcolorindex(L, 2, QC(V));
            mLatComplex *r = qlua_newLatComplex(L, Sidx);

            CALL_QDP(L);
            Qx(QDP_D,_C_eq_elem_V)(r->ptr, V->ptr, c, *S->qss);
        } else {
#if QNc == 'N'
            typedef QLA_DN_ColorVector(V->nc, Vtype);
#else
            typedef Qx(QLA_D,_ColorVector) Vtype;
#endif
            /* V[{x,y,z,...}], lattice coord must be present and valid */
            qlua_verifylatcoord(L, idx, S);
            int c = qlua_colorindex(L, 2, QC(V));
            CALL_QDP(L);
            Vtype *locked = Qx(QDP_D,_expose_V)(V->ptr);
            int site_node = QDP_node_number_L(S->lat, idx);

            if (c == -1) {
                /* V[{x,y,z,...}] -- get a color vector at a point */
                Qs(mSeqColVec) *v = Qs(qlua_newSeqColVec)(L, QC(V));
                double v_std[2 * QC(V)];
                
                if (site_node == QDP_this_node) {
                    int la = QDP_index_L(S->lat, idx);
                    Qs(unpack_vec)(QC(V), v_std, &locked[la]);
                }
                XMP_dist_double_array(site_node, 2 * QC(V), v_std);
                Qs(pack_vec)(QC(V), v->ptr, v_std);
            } else {
                /* V[{x,y,z,..., c = v}] -- one color component at a point */
                QLA_D_Complex *W = qlua_newComplex(L);
                double zri[2];
                qlua_checkcolorindex(L, 2, QC(V));

                if (site_node == QDP_this_node) {
                    int la = QDP_index_L(S->lat, idx);
                    QLA_D_Complex *z = &Qx(QLA_D,_elem_V)(locked[la], c);
                    zri[0] = QLA_real(*z);
                    zri[1] = QLA_imag(*z);
                }
                XMP_dist_double_array(site_node, 2, zri);
                QLA_c_eq_r_plus_ir(*W, zri[0], zri[1]);
            }
            Qx(QDP_D,_reset_V)(V->ptr);
            qlua_free(L, idx);
        }
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
    int *idx = qlua_latcoord(L, 2, S);

    if (idx == NULL) {
        /* V[{c=cx}] = cval */
        int c = qlua_checkcolorindex(L, 2, QC(V));
        mLatComplex *z = qlua_checkLatComplex(L, 3, S);
        
        CALL_QDP(L);
        Qx(QDP_D,_V_eq_elem_C)(V->ptr, z->ptr, c, *S->qss);
    } else {
#if QNc == 'N'
        typedef QLA_DN_ColorVector(V->nc, Vtype);
#else
        typedef Qx(QLA_D,_ColorVector) Vtype;
#endif
        int c = qlua_colorindex(L, 2, QC(V));
        qlua_verifylatcoord(L, idx, S);
        CALL_QDP(L);
        Vtype *locked = Qx(QDP_D,_expose_V)(V->ptr);
        int site_node = QDP_node_number_L(S->lat, idx);

        if (c == -1) {
            /* V[{x,y,...}] = v */
            Qs(mSeqColVec) *v = Qs(qlua_checkSeqColVec)(L, 3, QC(V));
            if (site_node == QDP_this_node) {
                int la = QDP_index_L(S->lat, idx);
                Qx(QLA_D,_V_eq_V)(QNC(QC(V)) &locked[la], v->ptr);
            }
        } else {
            /* V[{x,y,...,c=cx}] = c */
            qlua_checkcolorindex(L, 2, QC(V));
            QLA_D_Complex *z = qlua_checkComplex(L, 3);
            if (site_node == QDP_this_node) {
                int la = QDP_index_L(S->lat, idx);
                QLA_D_Complex *w = &Qx(QLA_D,_elem_V)(locked[la], c);
                QLA_c_eq_c(*w, *z);
            }
        }
        Qx(QDP_D,_reset_V)(V->ptr);
        qlua_free(L, idx);
    }
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

static int
Qs(q_R_mul_V_)(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColVec) *b = Qs(qlua_checkLatColVec)(L, 2, S, -1);
    Qs(mLatColVec) *r = Qs(qlua_newLatColVec)(L, lua_gettop(L), QC(b));

    CALL_QDP(L);
    Qx(QDP_D,_V_eq_R_times_V)(r->ptr, a->ptr, b->ptr, *S->qss);

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
    Qx(QDP_D,_V_eq_R_times_V)(r->ptr, b->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_C_mul_V_)(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColVec) *b = Qs(qlua_checkLatColVec)(L, 2, S, -1);
    Qs(mLatColVec) *r = Qs(qlua_newLatColVec)(L, lua_gettop(L), QC(b));

    CALL_QDP(L);
    Qx(QDP_D,_V_eq_C_times_V)(r->ptr, a->ptr, b->ptr, *S->qss);

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
    Qx(QDP_D,_V_eq_C_times_V)(r->ptr, b->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_V_sum)(lua_State *L)
{
    Qs(mLatColVec) *a = Qs(qlua_checkLatColVec)(L, 1, NULL, -1);
    int argc = lua_gettop(L);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    int nc = QC(a);

    switch (argc) {
    case 1: {
        Qs(mSeqColVec) *s = Qs(qlua_newSeqColVec)(L, nc);

        CALL_QDP(L);
        if (S->lss.mask) {
            Qs(mLatColVec) *b = Qs(qlua_newZeroLatColVec)(L, Sidx, nc);
            Qx(QDP_D,_V_eq_V_mask_I)(b->ptr, a->ptr, S->lss.mask, *S->qss);
            Qx(QDP_D,_v_eq_sum_V)(s->ptr, b->ptr, *S->qss);
            lua_pop(L, 1);
        } else {
            Qx(QDP_D,_v_eq_sum_V)(s->ptr, a->ptr, *S->qss);
        }
        return 1;
    }
    case 2: {
#if QNc == 'N'
        typedef QLA_DN_ColorVector(nc, Vtype);
#else
        typedef Qx(QLA_D,_ColorVector) Vtype;
#endif
        mLatMulti *m = qlua_checkLatMulti(L, 2, S);
        int size = m->size;
        QLA_Int *ii = m->idx;
        int sites = QDP_sites_on_node_L(S->lat);
        Vtype **vv = qlua_malloc(L, size * sizeof (Vtype *));
        int i, k, c;

        lua_createtable(L, size, 0);
        for (i = 0; i < size; i++) {
            Qs(mSeqColVec) *vi = Qs(qlua_newZeroSeqColVec)(L, QC(a));
            vv[i] = vi->ptr;
            lua_rawseti(L, -2, i + 1); /* [sic] lua index */
        }
        CALL_QDP(L);
        Vtype *xx = Qx(QDP_D,_expose_V)(a->ptr);
        
        for (k = 0; k < sites; k++, xx++, ii++) {
            int t = *ii;
            if ((t < 0) || (t >= size))
                continue;
            Qx(QLA_D,_V_peq_V)(QNC(nc) vv[t], xx);
        }
        Qx(QDP_D,_reset_V)(a->ptr);
        QLA_D_Real *rr = qlua_malloc(L, 2 * size * nc * sizeof (QLA_D_Real));
        for (i = 0; i < size; i++) {
            for (c = 0; c < nc; c++) {
                QLA_D_Complex *z = &Qx(QLA_D,_elem_V)(*vv[i], c);
                rr[2 * (c + i * nc)] = QLA_real(*z);
                rr[2 * (c + i * nc) + 1] = QLA_imag(*z);
            }
        }
        QMP_sum_double_array(rr, 2 * size * nc);
        for (i = 0; i < size; i++) {
            for (c = 0; c < nc; c++) {
                QLA_D_Complex z;
                QLA_c_eq_r_plus_ir(z, rr[2 *(c+i*nc)], rr[2*(c+i*nc)+1]);
                QLA_c_eq_c(Qx(QLA_D,_elem_V)(*vv[i], c), z);
            }
        }
                qlua_free(L, rr);
                qlua_free(L, vv);
        return 1;
    }
    }
    return luaL_error(L, "bad arguments for ColorVector:sum()");
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
    if (S->lss.mask) {
        Qx(QDP_D,_V_eq_V_mask_I)(r->ptr, a->ptr, S->lss.mask, *S->qss);
    } else {
        Qx(QDP_D,_V_eq_V)(r->ptr, a->ptr, *S->qss);
    }
    lua_pop(L, 2);

    return 0;
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
    { "sum",               Qs(q_V_sum)    },
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
    mLattice *S = qlua_checkLattice(L, Sidx);
#if QNc == 'N'
    Qx(QDP_D,_ColorVector) *v = Qx(QDP_D,_create_V_L)(nc, S->lat);
#else
    Qx(QDP_D,_ColorVector) *v = Qx(QDP_D,_create_V_L)(S->lat);
#endif
    Qs(mLatColVec) *hdr;

    if (v == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
#if QNc == 'N'
        v = Qx(QDP_D,_create_V_L)(nc, S->lat);
#else
        v = Qx(QDP_D,_create_V_L)(S->lat);
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
Qs(qlua_newZeroLatColVec)(lua_State *L, int Sidx, int nc)
{
        Qs(mLatColVec) *v = Qs(qlua_newLatColVec)(L, Sidx, nc);
        mLattice *S = qlua_checkLattice(L, Sidx);
        Qx(QDP_D,_V_eq_zero)(v->ptr, S->all);
        return v;
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
Qs(q_latcolvec_seq_)(lua_State *L, mLattice *S, int nc)
{
    Qs(mSeqColVec) *v = Qs(qlua_checkSeqColVec)(L, 2, nc);
    Qs(mLatColVec) *V = Qs(qlua_newLatColVec)(L, 1, nc);

    CALL_QDP(L);
    Qx(QDP_D, _V_eq_v)(V->ptr, v->ptr, *S->qss);
    
    return 1;
}

static int
Qs(q_latcolvec_)(lua_State *L, mLattice *S, int nc, int off)
{
    switch (lua_gettop(L) - off) {
    case 1: {
        Qs(qlua_newZeroLatColVec)(L, 1, nc);
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

static const QLUA_Op2 Qs(ops)[] = {
    { qlua_add_table, Qs(qLatColVec),  Qs(qLatColVec),  Qs(q_V_add_V_) },
    { qlua_sub_table, Qs(qLatColVec),  Qs(qLatColVec),  Qs(q_V_sub_V_) },
    { qlua_mul_table, qReal,           Qs(qLatColVec),  Qs(q_r_mul_V_) },
    { qlua_mul_table, Qs(qLatColVec),  qReal,           Qs(q_V_mul_r_) },
    { qlua_mul_table, qComplex,        Qs(qLatColVec),  Qs(q_c_mul_V_) },
    { qlua_mul_table, Qs(qLatColVec),  qComplex,        Qs(q_V_mul_c_) },
    { qlua_mul_table, qLatReal,        Qs(qLatColVec),  Qs(q_R_mul_V_) },
    { qlua_mul_table, Qs(qLatColVec),  qLatReal,        Qs(q_V_mul_R_) },
    { qlua_mul_table, qLatComplex,     Qs(qLatColVec),  Qs(q_C_mul_V_) },
    { qlua_mul_table, Qs(qLatColVec),  qLatComplex,     Qs(q_V_mul_C_) },
    { qlua_div_table, Qs(qLatColVec),  qReal,           Qs(q_V_div_r_) },
    { qlua_div_table, Qs(qLatColVec),  qComplex,        Qs(q_V_div_c_) },
    { NULL,           qNoType,         qNoType,         NULL           }
};

#undef QNc
#undef Qcolors
#undef Qs
#undef Qx
#undef QC
#undef QNC
