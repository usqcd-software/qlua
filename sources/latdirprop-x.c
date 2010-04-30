
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

    Qx(QDP_D,_destroy_P)(b->ptr);
    b->ptr = 0;

    return 0;
}

static void
#if QNc == 'N'
Qs(unpack_prop)(int nc, double *std, QLA_DN_DiracPropagator(nc, (*src)))
#else
Qs(unpack_prop)(int nc, double *std, Qx(QLA_D,_DiracPropagator) *src)
#endif
{
    for (int ci = 0; ci < nc; ci++) {
        for (int cj = 0; cj < nc; cj++) {
            for (int si = 0; si < QDP_Ns; si++) {
                for (int sj = 0; sj < QDP_Ns; sj++) {
                    QLA_D_Complex *z = &Qx(QLA_D,_elem_P)(*src, ci, si, cj, sj);
                    int la = 2 * (ci + nc * (cj + nc * (si + QDP_Ns * sj)));
                    std[la] = QLA_real(*z);
                    std[la + 1] = QLA_imag(*z);
                }
            }
        }
    }
}

static void
#if QNc == 'N'
Qs(pack_prop)(int nc, QLA_DN_DiracPropagator(nc, (*dst)), const double *std)
#else
Qs(pack_prop)(int nc, Qx(QLA_D,_DiracPropagator) *dst, const double *std)
#endif
{
    QLA_D_Complex z;
    for (int ci = 0; ci < nc; ci++) {
        for (int cj = 0; cj < nc; cj++) {
            for (int si = 0; si < QDP_Ns; si++) {
                for (int sj = 0; sj < QDP_Ns; sj++) {
                    int la = 2 * (ci + nc * (cj + nc * (si + QDP_Ns * sj)));
                    QLA_c_eq_r_plus_ir(z, std[la], std[la + 1]);
                    QLA_c_eq_c(Qx(QLA_D,_elem_P)(*dst, ci, si, cj, sj), z);
                }
            }
        }
    }
}

static int
Qs(q_P_get)(lua_State *L)
{
    switch (qlua_qtype(L, 2)) {
    case qTable: {
        Qs(mLatDirProp) *V = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
        mLattice *S = qlua_ObjLattice(L, 1);
        int Sidx = lua_gettop(L);
        int d = qlua_diracindex(L, 2);
        int c = qlua_colorindex(L, 2, QC(V));
        int *idx = qlua_latcoord(L, 2, S);
#if QNc == 'N'
        typedef QLA_DN_DiracPropagator(QC(V), Vtype);
#else
        typedef Qx(QLA_D,_DiracPropagator) Vtype;
#endif

        if (idx == 0) {
            /* P[{d=expr,c=expr}] => D */
            if ((d == -1) || (c == -1))
                return qlua_badindex(L, "DiracPropagator" Qcolors);
            Qs(mLatDirFerm) *r = Qs(qlua_newLatDirFerm)(L, Sidx, QC(V));
                
            CALL_QDP(L);
            Qx(QDP_D,_D_eq_diracvec_P)(r->ptr, V->ptr, c, d, *S->qss);
        } else {
            /* P[{x,...}] => p */
            qlua_verifylatcoord(L, idx, S);
            if ((d != -1) || (c != -1))
                return qlua_badindex(L, "DiracPropagator" Qcolors);
            Qs(mSeqDirProp) *m = Qs(qlua_newSeqDirProp)(L, QC(V));
            int len = 2 * QC(V) * QC(V) * QDP_Ns * QDP_Ns;
            double m_std[len];
            CALL_QDP(L);
            Vtype *locked = Qx(QDP_D,_expose_P)(V->ptr);
            int site_node = QDP_node_number_L(S->lat, idx);
            if (site_node == QDP_this_node) {
                int la = QDP_index_L(S->lat, idx);
                Qs(unpack_prop)(QC(V), m_std, &locked[la]);
            }
            XMP_dist_double_array(site_node, len, m_std);
            Qs(pack_prop)(QC(V), m->ptr, m_std);
            Qx(QDP_D,_reset_P)(V->ptr);
            qlua_free(L, idx);
        }
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
    int d = qlua_diracindex(L, 2);
    int c = qlua_colorindex(L, 2, QC(V));
    int *idx = qlua_latcoord(L, 2, S);

    if (idx == NULL) {
        /* P[{c=expr,d=expr}] = D */
        if ((d == -1) || (c == -1))
            return qlua_badindex(L, "DiracPropagator" Qcolors);
        Qs(mLatDirFerm) *r = Qs(qlua_checkLatDirFerm)(L, 3, S, QC(V));
            
        CALL_QDP(L);
        Qx(QDP_D,_P_eq_diracvec_D)(V->ptr, r->ptr, c, d, *S->qss);
    } else {
        /* P[{x,y,z...}] = p */
        qlua_verifylatcoord(L, idx, S);
        if ((d != -1) || (c != -1))
            return qlua_badindex(L, "DiracPropagator" Qcolors);
#if QNc == 'N'
        typedef QLA_DN_DiracPropagator(QC(V), Vtype);
#else
        typedef Qx(QLA_D,_DiracPropagator) Vtype;
#endif
        CALL_QDP(L);
        Vtype *locked = Qx(QDP_D,_expose_P)(V->ptr);
        int site_node = QDP_node_number_L(S->lat, idx);
        Qs(mSeqDirProp) *v = Qs(qlua_checkSeqDirProp)(L, 3, QC(V));
        if (site_node == QDP_this_node) {
            int la = QDP_index_L(S->lat, idx);
            Qx(QLA_D,_P_eq_P)(QNC(QC(V)) &locked[la], v->ptr);
        }
        Qx(QDP_D,_reset_P)(V->ptr);
        qlua_free(L, idx);
    }
    return 0;
}

static int
Qs(q_P_sum)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    int argc = lua_gettop(L);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    int nc = QC(a);
    
    switch (argc) {
    case 1: {
        Qs(mSeqDirProp) *s = Qs(qlua_newSeqDirProp)(L, nc);

        CALL_QDP(L);
        if (S->lss.mask) {
            Qs(mLatDirProp) *b = Qs(qlua_newLatDirProp)(L, Sidx, nc);
            Qx(QDP_D,_P_eq_zero)(b->ptr, *S->qss);
            Qx(QDP_D,_P_eq_P_mask_I)(b->ptr, a->ptr, S->lss.mask, *S->qss);
            Qx(QDP_D,_p_eq_sum_P)(s->ptr, b->ptr, *S->qss);
            lua_pop(L, 1);
        } else {
            Qx(QDP_D,_p_eq_sum_P)(s->ptr, a->ptr, *S->qss);
        }
        return 1;
    }
    case 2: {
#if QNc == 'N'
        typedef QLA_DN_DiracPropagator(nc, Vtype);
#else
        typedef Qx(QLA_D,_DiracPropagator) Vtype;
#endif
        mLatMulti *m = qlua_checkLatMulti(L, 2, S);
        int size = m->size;
        QLA_Int *ii = m->idx;
        int sites = QDP_sites_on_node_L(S->lat);

        Vtype *vv[size];
        lua_createtable(L, size, 0);
        for (int i = 0; i < size; i++) {
            Qs(mSeqDirProp) *vi = Qs(qlua_newSeqDirProp)(L, QC(a));
            Qx(QLA_D,_P_eq_zero)(QNC(nc) vi->ptr);
            vv[i] = vi->ptr;
            lua_rawseti(L, -2, i + 1); /* [sic] lua index */
        }
        CALL_QDP(L);
        Vtype *xx = Qx(QDP_D,_expose_P)(a->ptr);
        
        for (int k = 0; k < sites; k++, xx++, ii++) {
            int t = *ii;
            if ((t < 0) || (t >= size))
                continue;
            Qx(QLA_D,_P_peq_P)(QNC(nc) vv[t], xx);
        }
        Qx(QDP_D,_reset_P)(a->ptr);
        QLA_D_Real rr[2 * size * nc * nc * QDP_Ns * QDP_Ns];
        for (int i = 0; i < size; i++) {
            for (int ci = 0; ci < nc; ci++) {
                for (int cj = 0; cj < nc; cj++) {
                    for (int di = 0; di < QDP_Ns; di++) {
                        for (int dj = 0; dj < QDP_Ns; dj++) {
                            QLA_D_Complex *z;
                            int ab = 2*(ci + nc*(cj + nc*(di + QDP_Ns*dj)));
                            z = &Qx(QLA_D,_elem_P)(*vv[i], ci, di, cj, dj);
                            rr[ab] = QLA_real(*z);
                            rr[ab + 1] = QLA_imag(*z);
                        }
                    }
                }
            }
        }
        QMP_sum_double_array(rr, 2 * size * nc * nc * QDP_Ns * QDP_Ns);
        for (int i = 0; i < size; i++) {
            for (int ci = 0; ci < nc; ci++) {
                for (int cj = 0; cj < nc; cj++) {
                    for (int di = 0; di < QDP_Ns; di++) {
                        for (int dj = 0; dj < QDP_Ns; dj++) {
                            QLA_D_Complex z;
                            int ab = 2*(ci + nc*(cj + nc*(di + QDP_Ns*dj)));
                            QLA_c_eq_r_plus_ir(z, rr[ab], rr[ab + 1]);
                            QLA_c_eq_c(Qx(QLA_D,_elem_P)(*vv[i],ci,di,cj,dj),z);
                        }
                    }
                }
            }
        }
        return 1;
    }
    }
    return luaL_error(L, "bad arguments for DiracPropagator:sum()");
}

static int
Qs(q_P_norm2_)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    QLA_D_Real n;

    CALL_QDP(L);
    if (S->lss.mask) {
        Qs(mLatDirProp) *b = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(a));
        Qx(QDP_D,_P_eq_P_mask_I)(b->ptr, a->ptr, S->lss.mask, *S->qss);
        Qx(QDP_D,_r_eq_norm2_P)(&n, b->ptr, *S->qss);
    } else {
        Qx(QDP_D,_r_eq_norm2_P)(&n, a->ptr, *S->qss);
    }
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
    Qx(QDP_D,_P_eq_sP)(r->ptr, a->ptr, shift, dir, *S->qss);

    return 1;
}

static int
Qs(q_P_conj)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirProp) *r = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_P_eq_conj_P)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_P_trans)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirProp) *r = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_P_eq_transpose_P)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_P_adjoin)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirProp) *r = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_P_eq_Pa)(r->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_P_spintrace)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatColMat) *r = Qs(qlua_newLatColMat)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_M_eq_spintrace_P)(r->ptr, a->ptr, *S->qss);

    return 1;
}

typedef struct {
    int nc;
    void *a; /* Qx(QLA_D,_DiracPropagator) *a; */
} Qs(Pst_args);

static void
#if QNc == 'N'
Qs(do_Pst)(int nc, QLA_DN_DiracPropagator(nc, (*r)), int idx, void *env)
#else
Qs(do_Pst)(Qx(QLA_D,_DiracPropagator) *r, int idx, void *env)
#endif
{
    Qs(Pst_args) *arg = env;
#if QNc == 'N'
    QLA_DN_DiracPropagator(nc, (*ai)) = arg->a;
    QLA_DN_DiracPropagator(nc, (*a)) = &ai[idx];
#else
    Qx(QLA_D,_DiracPropagator) *ai = arg->a;
    Qx(QLA_D,_DiracPropagator) *a = &ai[idx];
    int nc = QC(a);
#endif
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
    Qs(Pst_args) arg;

    CALL_QDP(L);
    arg.nc = QC(a);
    arg.a = Qx(QDP_D,_expose_P)(a->ptr);
    Qx(QDP_D,_P_eq_funcia)(r->ptr, Qs(do_Pst), &arg, *S->qss);
    Qx(QDP_D,_reset_P)(a->ptr);

    return 1;
}

static void
Qs(do_Ptrace)(QLA_D_Complex *r, int idx, void *env)
{
    Qs(Pst_args) *arg = env;
    int nc = arg->nc;
#if QNc == 'N'
    QLA_DN_DiracPropagator(nc, (*ai)) = arg->a;
    QLA_DN_DiracPropagator(nc, (*a)) = &ai[idx];
#else
    Qx(QLA_D,_DiracPropagator) *ai = arg->a;
    Qx(QLA_D,_DiracPropagator) *a = &ai[idx];
#endif
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
    Qs(Pst_args) arg;

    CALL_QDP(L);
    arg.nc = QC(a);
    arg.a = Qx(QDP_D,_expose_P)(a->ptr);
    QDP_D_C_eq_funcia(r->ptr, Qs(do_Ptrace), &arg, *S->qss);
    Qx(QDP_D,_reset_P)(a->ptr);

    return 1;
}

static int
Qs(q_P_set)(lua_State *L)
{
    Qs(mLatDirProp) *r = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 2, S, QC(r));

    CALL_QDP(L);
    if (S->lss.mask)
        Qx(QDP_D,_P_eq_P_mask_I)(r->ptr, a->ptr, S->lss.mask, *S->qss);
    else
        Qx(QDP_D,_P_eq_P)(r->ptr, a->ptr, *S->qss);

    return 0;
}

static int
Qs(q_P_add_P_)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirProp) *b = Qs(qlua_checkLatDirProp)(L, 2, S, QC(a));
    Qs(mLatDirProp) *c = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_P_eq_P_plus_P)(c->ptr, a->ptr, b->ptr, *S->qss);

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
    Qx(QDP_D,_P_eq_P_minus_P)(c->ptr, a->ptr, b->ptr, *S->qss);

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
    Qx(QDP_D,_P_eq_r_times_P)(c->ptr, &b, a->ptr, *S->qss);

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
    Qx(QDP_D,_P_eq_r_times_P)(c->ptr, &a, b->ptr, *S->qss);

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
    Qx(QDP_D,_P_eq_c_times_P)(c->ptr, b, a->ptr, *S->qss);

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
    Qx(QDP_D,_P_eq_c_times_P)(c->ptr, a, b->ptr, *S->qss);

    return 1;
}

static int
Qs(q_R_mul_P_)(lua_State *L)
{
    mLatReal *a = qlua_checkLatReal(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirProp) *b = Qs(qlua_checkLatDirProp)(L, 2, S, -1);
    Qs(mLatDirProp) *r = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(b));

    CALL_QDP(L);
    Qx(QDP_D,_P_eq_R_times_P)(r->ptr, a->ptr, b->ptr, *S->qss);

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
    Qx(QDP_D,_P_eq_R_times_P)(r->ptr, b->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_C_mul_P_)(lua_State *L)
{
    mLatComplex *a = qlua_checkLatComplex(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirProp) *b = Qs(qlua_checkLatDirProp)(L, 2, S, -1);
    Qs(mLatDirProp) *r = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(b));

    CALL_QDP(L);
    Qx(QDP_D,_P_eq_C_times_P)(r->ptr, a->ptr, b->ptr, *S->qss);

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
    Qx(QDP_D,_P_eq_C_times_P)(r->ptr, b->ptr, a->ptr, *S->qss);

    return 1;
}

static int
Qs(q_P_neg)(lua_State *L)
{
    Qs(mLatDirProp) *a = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    Qs(mLatDirProp) *r = Qs(qlua_newLatDirProp)(L, lua_gettop(L), QC(a));

    CALL_QDP(L);
    Qx(QDP_D,_P_meq_P)(r->ptr, a->ptr, *S->qss);

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
    Qx(QDP_D,_P_eq_P_times_P)(c->ptr, a->ptr, b->ptr, *S->qss);

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
    Qx(QDP_D,_P_eq_P_times_M)(c->ptr, a->ptr, b->ptr, *S->qss);

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
    Qx(QDP_D,_P_eq_M_times_P)(c->ptr, a->ptr, b->ptr, *S->qss);

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
    Qx(QDP_D,_P_eq_r_times_P)(c->ptr, &b, a->ptr, *S->qss);

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
    Qx(QDP_D,_P_eq_c_times_P)(c->ptr, &s, a->ptr, *S->qss);

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
    Qx(QDP_D,_C_eq_P_dot_P)(s->ptr, a->ptr, b->ptr, *S->qss);

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
    Qx(QDP_D,_P_eq_P)(r->ptr, a->ptr, *S->qss);

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
    { "sum",               Qs(q_P_sum)            },
    { "norm2",             Qs(q_P_norm2_)         },
    { "shift",             Qs(q_P_shift)          },
    { "conj",              Qs(q_P_conj)           },
    { "transpose",         Qs(q_P_trans)          },
    { "adjoin",            Qs(q_P_adjoin)         },
    { "spintrace",         Qs(q_P_spintrace)      },
    { "spintranspose",     Qs(q_P_spintranspose)  },
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
    mLattice *S = qlua_checkLattice(L, Sidx);
#if QNc == 'N'
    Qx(QDP_D,_DiracPropagator) *v = Qx(QDP_D,_create_P_L)(nc, S->lat);
#else
    Qx(QDP_D,_DiracPropagator) *v = Qx(QDP_D,_create_P_L)(S->lat);
#endif
    Qs(mLatDirProp) *hdr;

    if (v == 0) {
        lua_gc(L, LUA_GCCOLLECT, 0);
#if QNc == 'N'
        v = Qx(QDP_D,_create_P_L)(nc, S->lat);
#else
        v = Qx(QDP_D,_create_P_L)(S->lat);
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
Qs(q_latdirprop_seq_)(lua_State *L, mLattice *S, int nc)
{
    Qs(mSeqDirProp) *m = Qs(qlua_checkSeqDirProp)(L, 2, nc);
    Qs(mLatDirProp) *M = Qs(qlua_newLatDirProp)(L, 1, nc);

    CALL_QDP(L);
    Qx(QDP_D, _P_eq_p)(M->ptr, m->ptr, *S->qss);

    return 1;
}

static int
Qs(q_latdirprop_lat_)(lua_State *L, mLattice *S, int nc)
{
    Qs(mLatDirProp) *m = Qs(qlua_checkLatDirProp)(L, 2, S, nc);
    Qs(mLatDirProp) *M = Qs(qlua_newLatDirProp)(L, 1, nc);

    CALL_QDP(L);
    Qx(QDP_D, _P_eq_P)(M->ptr, m->ptr, *S->qss);

    return 1;
}

static int
Qs(q_latdirprop_mat_)(lua_State *L, mLattice *S, int nc)
{
    Qs(mLatColMat) *w = Qs(qlua_checkLatColMat)(L, 2, S, nc);
    Qs(mLatDirProp) *v = Qs(qlua_newLatDirProp)(L, 1, nc);
    mLatComplex *c = qlua_newLatComplex(L, 1);
    int ic, jc, ks;

    CALL_QDP(L);
    Qx(QDP_D,_P_eq_zero)(v->ptr, *S->qss);
    for (ic = 0; ic < nc; ic++) {
        for (jc = 0; jc < nc; jc++) {
            Qx(QDP_D,_C_eq_elem_M)(c->ptr, w->ptr, ic, jc, *S->qss);
            for (ks = 0; ks < QDP_Ns; ks++)
                Qx(QDP_D,_P_eq_elem_C)(v->ptr, c->ptr,
                                       ic, ks, jc, ks, *S->qss);
        }
    }
    lua_pop(L, 1);
    return 1;
}

static int
Qs(q_latdirprop_)(lua_State *L, mLattice *S, int nc, int off)
{
    switch (lua_gettop(L) - off) {
    case 1: {
        Qs(mLatDirProp) *v = Qs(qlua_newLatDirProp)(L, 1, nc);

        CALL_QDP(L);
        Qx(QDP_D,_P_eq_zero)(v->ptr, *S->qss);
        
        return 1;
    }
    case 3: {
        Qs(mLatDirFerm) *z = Qs(qlua_checkLatDirFerm)(L, 2 + off, S, nc);
        int d = qlua_checkdiracindex(L, 3);
        int c = qlua_checkcolorindex(L, 3, nc);
        Qs(mLatDirProp) *v = Qs(qlua_newLatDirProp)(L, 1, nc);

        CALL_QDP(L);
        Qx(QDP_D,_P_eq_zero)(v->ptr, *S->qss);
        Qx(QDP_D,_P_eq_diracvec_D)(v->ptr, z->ptr, c, d, *S->qss);
        return 1;
    }
    }
    return qlua_badconstr(L, "DiracPropagator" Qcolors);
}
static const QLUA_Op2 Qs(ops)[] = {
    { qlua_add_table, Qs(qLatDirProp),  Qs(qLatDirProp),  Qs(q_P_add_P_) },
    { qlua_sub_table, Qs(qLatDirProp),  Qs(qLatDirProp),  Qs(q_P_sub_P_) },
    { qlua_mul_table, qReal,            Qs(qLatDirProp),  Qs(q_r_mul_P_) },
    { qlua_mul_table, Qs(qLatDirProp),  qReal,            Qs(q_P_mul_r_) },
    { qlua_mul_table, qComplex,         Qs(qLatDirProp),  Qs(q_c_mul_P_) },
    { qlua_mul_table, Qs(qLatDirProp),  qComplex,         Qs(q_P_mul_c_) },
    { qlua_mul_table, qLatReal,         Qs(qLatDirProp),  Qs(q_R_mul_P_) },
    { qlua_mul_table, Qs(qLatDirProp),  qLatReal,         Qs(q_P_mul_R_) },
    { qlua_mul_table, qLatComplex,      Qs(qLatDirProp),  Qs(q_C_mul_P_) },
    { qlua_mul_table, Qs(qLatDirProp),  qLatComplex,      Qs(q_P_mul_C_) },
    { qlua_mul_table, Qs(qLatDirProp),  Qs(qLatDirProp),  Qs(q_P_mul_P_) },
    { qlua_mul_table, Qs(qLatDirProp),  Qs(qLatColMat),   Qs(q_P_mul_M_) },
    { qlua_mul_table, Qs(qLatColMat),   Qs(qLatDirProp),  Qs(q_M_mul_P_) },
    { qlua_div_table, Qs(qLatDirProp),  qReal,            Qs(q_P_div_r_) },
    { qlua_div_table, Qs(qLatDirProp),  qComplex,         Qs(q_P_div_c_) },
    { NULL,           qNoType,          qNoType,          NULL           }
};

#undef QNc
#undef Qcolors
#undef Qs
#undef Qx
#undef QC
#undef QNC
