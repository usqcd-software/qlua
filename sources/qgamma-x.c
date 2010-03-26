#define PIDX(a,b)  ((a) * QDP_Ns/2 + (b))
static int
Qs(q_g_mul_D_)(lua_State *L)
{
    mClifford *m = qlua_checkClifford(L, 1);
    Qs(mLatDirFerm) *f = Qs(qlua_checkLatDirFerm)(L, 2, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 2);
    int Sidx = lua_gettop(L);
    Qs(mLatDirFerm) *mf = Qs(qlua_newLatDirFerm)(L, Sidx, QC(f));
    Qs(mLatDirFerm) *r = Qs(qlua_newLatDirFerm)(L, Sidx, QC(f));
    int i;

    CALL_QDP(L);
    Qx(QDP_D,_D_eq_zero)(r->ptr, *S->qss);
    for (i = 0; i < 16; i++) {
        switch (m->g[i].t) {
        case qG_z: continue;
        case qG_p:
            Qx(QDP_D,_D_eq_gamma_times_D)(mf->ptr, f->ptr, i, *S->qss);
            Qx(QDP_D,_D_peq_D)(r->ptr, mf->ptr, *S->qss);
            break;
        case qG_m:
            Qx(QDP_D,_D_eq_gamma_times_D)(mf->ptr, f->ptr, i, *S->qss);
            Qx(QDP_D,_D_meq_D)(r->ptr, mf->ptr, *S->qss);
            break;
        case qG_r:
            Qx(QDP_D,_D_eq_gamma_times_D)(mf->ptr, f->ptr, i, *S->qss);
            Qx(QDP_D,_D_peq_r_times_D)(r->ptr, &m->g[i].r, mf->ptr, *S->qss);
            break;
        case qG_c:
            Qx(QDP_D,_D_eq_gamma_times_D)(mf->ptr, f->ptr, i, *S->qss);
            Qx(QDP_D,_D_peq_c_times_D)(r->ptr, &m->g[i].c, mf->ptr, *S->qss);
            break;
        }
    }
    return 1;
}

static int
Qs(q_g_mul_P_)(lua_State *L)
{
    mClifford *m = qlua_checkClifford(L, 1);
    Qs(mLatDirProp) *f = Qs(qlua_checkLatDirProp)(L, 2, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 2);
    int Sidx = lua_gettop(L);
    Qs(mLatDirProp) *mf = Qs(qlua_newLatDirProp)(L, Sidx, QC(f));
    Qs(mLatDirProp) *r = Qs(qlua_newLatDirProp)(L, Sidx, QC(f));
    int i;
    
    CALL_QDP(L);
    Qx(QDP_D,_P_eq_zero)(r->ptr, *S->qss);
    for (i = 0; i < 16; i++) {
        switch (m->g[i].t) {
        case qG_z: continue;
        case qG_p:
            Qx(QDP_D,_P_eq_gamma_times_P)(mf->ptr, f->ptr, i, *S->qss);
            Qx(QDP_D,_P_peq_P)(r->ptr, mf->ptr, *S->qss);
            break;
        case qG_m:
            Qx(QDP_D,_P_eq_gamma_times_P)(mf->ptr, f->ptr, i, *S->qss);
            Qx(QDP_D,_P_meq_P)(r->ptr, mf->ptr, *S->qss);
            break;
        case qG_r:
            Qx(QDP_D,_P_eq_gamma_times_P)(mf->ptr, f->ptr, i, *S->qss);
            Qx(QDP_D,_P_peq_r_times_P)(r->ptr, &m->g[i].r, mf->ptr, *S->qss);
            break;
        case qG_c:
               Qx(QDP_D,_P_eq_gamma_times_P)(mf->ptr, f->ptr, i, *S->qss);
               Qx(QDP_D,_P_peq_c_times_P)(r->ptr, &m->g[i].c, mf->ptr, *S->qss);
            break;
        }
    }
    return 1;
}

static int
Qs(q_P_mul_g_)(lua_State *L)
{
    Qs(mLatDirProp) *f = Qs(qlua_checkLatDirProp)(L, 1, NULL, -1);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
    mClifford *m = qlua_checkClifford(L, 2);
    Qs(mLatDirProp) *mf = Qs(qlua_newLatDirProp)(L, Sidx, QC(f));
    Qs(mLatDirProp) *r = Qs(qlua_newLatDirProp)(L, Sidx, QC(f));
    int i;

    CALL_QDP(L);
    Qx(QDP_D,_P_eq_zero)(r->ptr, *S->qss);
    for (i = 0; i < 16; i++) {
        switch (m->g[i].t) {
        case qG_z: continue;
        case qG_p:
            Qx(QDP_D,_P_eq_P_times_gamma)(mf->ptr, f->ptr, i, *S->qss);
            Qx(QDP_D,_P_peq_P)(r->ptr, mf->ptr, *S->qss);
            break;
        case qG_m:
            Qx(QDP_D,_P_eq_P_times_gamma)(mf->ptr, f->ptr, i, *S->qss);
            Qx(QDP_D,_P_meq_P)(r->ptr, mf->ptr, *S->qss);
            break;
        case qG_r:
            Qx(QDP_D,_P_eq_P_times_gamma)(mf->ptr, f->ptr, i, *S->qss);
            Qx(QDP_D,_P_peq_r_times_P)(r->ptr, &m->g[i].r, mf->ptr, *S->qss);
            break;
        case qG_c:
            Qx(QDP_D,_P_eq_P_times_gamma)(mf->ptr, f->ptr, i, *S->qss);
            Qx(QDP_D,_P_peq_c_times_P)(r->ptr, &m->g[i].c, mf->ptr, *S->qss);
            break;
        }
    }
    return 1;
}

#if QNc == 'N'
typedef void Qs(X_project)(int nc,
                  lua_State *L,
                  QLA_DN_DiracFermion(nc, (**r)),
                  int mu,
                  int sign,
                  QLA_DN_DiracPropagator(nc, (*f)));
#else
typedef void Qs(X_project)(int nc,
                  lua_State *L,
                  Qx(QLA_D,_DiracFermion) **r,
                  int mu,
                  int sign,
                  Qx(QLA_D,_DiracPropagator) *f);
#endif

static void
Qs(q_P_left_proj)(int nc,
                  lua_State *L,
#if QNc == 'N'                  
                  QLA_DN_DiracFermion(nc, (**r)),
#else
                  Qx(QLA_D,_DiracFermion) **r,
#endif            
                  int mu,
                  int sign,
#if QNc == 'N'
                  QLA_DN_DiracPropagator(nc, (*f))
#else
                  Qx(QLA_D,_DiracPropagator) *f
#endif
    )
{
#if QNc == 'N'
    typedef QLA_DN_DiracFermion(nc, Vtype);
    typedef QLA_DN_HalfFermion(nc, Htype);
#else
    typedef Qx(QLA_D,_DiracFermion) Vtype;
    typedef Qx(QLA_D,_HalfFermion) Htype;
#endif
    int count = QDP_sites_on_node;
    int k, ic, is;

    for (k = 0; k < count; k++) {
        for (ic = 0; ic < nc; ic++) {
            for (is = 0; is < QDP_Ns; is++) {
                int jc, js;
                Vtype fk;
                Htype hk;
                for (jc = 0; jc < nc; jc++) {
                    for (js = 0; js < QDP_Ns; js++) {
                        QLA_c_eq_c(QLA_elem_D(fk, jc, js),
                                   QLA_elem_P(f[0], jc, js, ic, is));
                    }
                }
#if QNc == 'N'
                Qx(QLA_D,_H_eq_spproj_D)(nc, &hk, &fk, mu, sign);
#else
                Qx(QLA_D,_H_eq_spproj_D)(&hk, &fk, mu, sign);
#endif          
                for (jc = 0; jc < nc; jc++) {
                    for (js = 0; js < QDP_Ns/2; js++) {
                        QLA_c_eq_c(QLA_elem_D(*r[PIDX(jc, js)], ic, is),
                                   QLA_elem_H(hk, jc, js));
                    }
                }
            }
        }
        f++;
        for (ic = 0; ic < nc; ic++)
            for (is = 0; is < QDP_Ns/2; is++)
                r[PIDX(ic, is)]++;
    }
}

static void
Qs(q_P_right_proj)(int nc,
                   lua_State *L,
#if QNc == 'N'
                   QLA_DN_DiracFermion(nc, (**r)),
#else
                   Qx(QLA_D,_DiracFermion) **r,
#endif
                   int mu,
                   int sign,
#if QNc == 'N'
                   QLA_DN_DiracPropagator(nc, (*f))
#else
                   Qx(QLA_D,_DiracPropagator) *f
#endif
    )
{
#if QNc == 'N'
    typedef QLA_DN_DiracFermion(nc, Vtype);
    typedef QLA_DN_HalfFermion(nc, Htype);
#else
    typedef Qx(QLA_D,_DiracFermion) Vtype;
    typedef Qx(QLA_D,_HalfFermion) Htype;
#endif
    int count = QDP_sites_on_node;
    int k, ic, is;

    if (mu == 0 || mu == 2) sign = 1 - sign; /* AAA gamma basis dependent */

    for (k = 0; k < count; k++) {
        for (ic = 0; ic < nc; ic++) {
            for (is = 0; is < QDP_Ns; is++) {
                int jc, js;
                Vtype fk;
                Htype hk;
                for (jc = 0; jc < nc; jc++) {
                    for (js = 0; js < QDP_Ns; js++) {
                        QLA_c_eq_c(QLA_elem_D(fk, jc, js),
                                   QLA_elem_P(*f, ic, is, jc, js));
                    }
                }
#if QNc == 'N'
                Qx(QLA_D,_H_eq_spproj_D)(nc, &hk, &fk, mu, sign);
#else
                Qx(QLA_D,_H_eq_spproj_D)(&hk, &fk, mu, sign);
#endif
                for (jc = 0; jc < nc; jc++) {
                    for (js = 0; js < QDP_Ns/2; js++) {
                        QLA_c_eq_c(QLA_elem_D(*r[PIDX(jc, js)], ic, is),
                                   QLA_elem_H(hk, jc, js));
                    }
                }
            }
        }
        f++;
        for (ic = 0; ic < nc; ic++)
            for (is = 0; is < QDP_Ns/2; is++)
                r[PIDX(ic, is)]++;
    }
}

static int
Qs(q_P_project)(lua_State *L, Qs(X_project) *op)
{
    const char *sign = luaL_checkstring(L, 1);
    int mu = qlua_checkgammaindex(L, 2);
    Qs(mLatDirProp) *f = Qs(qlua_checkLatDirProp)(L, 3, NULL, -1);
    int Sidx;
    int isign = 0;
    int ic, is;
    Qs(mLatDirFerm) *r[QC(f) * QDP_Ns / 2];
#if QNc == 'N'
    QLA_DN_DiracPropagator(f->nc, (*ff));
    QLA_DN_DiracFermion(f->nc, (*rr[QC(f) * QDP_Ns/2]));
#else
    Qx(QLA_D,_DiracPropagator) *ff;
    Qx(QLA_D,_DiracFermion) *rr[QC(f) * QDP_Ns/2];
#endif

    qlua_ObjLattice(L, 3);
    Sidx = lua_gettop(L);

    if (strcmp(sign, "plus") == 0)
        isign = 1;
    else if (strcmp(sign, "minus") == 0)
        isign = 0;
    else
        return luaL_error(L, "bad sign parameter");
    
    lua_createtable(L, QC(f), 0);
    for (ic = 0; ic < QC(f); ic++) {
        lua_createtable(L, QDP_Ns/2, 0);
        for (is = 0; is < QDP_Ns/2; is++) {
            r[PIDX(ic, is)] = Qs(qlua_newLatDirFerm)(L, Sidx, QC(f));
            lua_rawseti(L, -2, is + 1);
        }
        lua_rawseti(L, -2, ic + 1);
    }

    if (mu == 5) mu = 4; /* [sic] -- funny numbering in QLA's spproj/sprecon */

    CALL_QDP(L);
    for (ic = 0; ic < QC(f); ic++) {
        for (is = 0; is < QDP_Ns/2; is++) {
            int idx_cs = ic * QDP_Ns/2 + is;
            rr[idx_cs] = Qx(QDP_D,_expose_D)(r[idx_cs]->ptr);
        }
    }
    ff = Qx(QDP_D,_expose_P)(f->ptr);

    op(QC(f), L, rr, mu, isign, ff);

    Qx(QDP_D,_reset_P)(f->ptr);
    for (ic = 0; ic < QC(f); ic++) {
        for (is = 0; is < QDP_Ns/2; is++)  {
            int idx_cs = PIDX(ic, is);
            Qx(QDP_D,_reset_D)(r[idx_cs]->ptr);
        }
    }

    return 1;
}

#if QNc == 'N'
typedef void Qs(X_reconstruct)(int nc,
                               lua_State *L,
                               QLA_DN_DiracPropagator(nc, (*r)),
                               int mu, int sign,
                               QLA_DN_DiracFermion(nc, (**a)));
#else
typedef void Qs(X_reconstruct)(int nc,
                               lua_State *L,
                               Qx(QLA_D,_DiracPropagator) *r,
                               int mu, int sign,
                               Qx(QLA_D,_DiracFermion) **a);
#endif

#if QNc == 'N'
static void
Qs(q_P_left_recon)(int nc,
                   lua_State *L,
                   QLA_DN_DiracPropagator(nc, (*r)),
                   int mu, int sign,
                   QLA_DN_DiracFermion(nc, (**a)))
#else
static void
Qs(q_P_left_recon)(int nc,
                   lua_State *L,
                   Qx(QLA_D,_DiracPropagator) *r,
                   int mu, int sign,
                   Qx(QLA_D,_DiracFermion) **a)
#endif
{
#if QNc == 'N'
    typedef QLA_DN_DiracFermion(nc, Vtype);
    typedef QLA_DN_HalfFermion(nc, Htype);
#else
    typedef Qx(QLA_D,_DiracFermion) Vtype;
    typedef Qx(QLA_D,_HalfFermion) Htype;
#endif
    int count = QDP_sites_on_node;
    int k, ic, is;

    for (k = 0; k < count; k++) {
        for (ic = 0; ic < nc; ic++) {
            for (is = 0; is < QDP_Ns; is++) {
                int jc, js;
                Vtype fk;
                Htype hk;
                for (jc = 0; jc < nc; jc++) {
                    for (js = 0; js < QDP_Ns/2; js++) {
                        QLA_c_eq_c(QLA_elem_H(hk, jc, js),
                                   QLA_elem_D(*a[PIDX(jc,js)], ic, is));
                    }
                }
#if QNc == 'N'
                Qx(QLA_D,_D_eq_sprecon_H)(nc, &fk, &hk, mu, sign);
#else
                Qx(QLA_D,_D_eq_sprecon_H)(&fk, &hk, mu, sign);
#endif
                for (jc = 0; jc < nc; jc++) {
                    for (js = 0; js < QDP_Ns; js++) {
                        QLA_c_eq_c(QLA_elem_P(*r, jc, js, ic, is),
                                   QLA_elem_D(fk, jc, js));
                    }
                }
            }
        }
        r++;
        for (ic = 0; ic < nc; ic++)
            for (is = 0; is < QDP_Ns/2; is++)
                a[PIDX(ic, is)]++;
    }
    
}

#if QNc == 'N'
static void
Qs(q_P_right_recon)(int nc,
                    lua_State *L,
                    QLA_DN_DiracPropagator(nc, (*r)),
                    int mu, int sign,
                    QLA_DN_DiracFermion(nc, (**a)))
#else
static void
Qs(q_P_right_recon)(int nc,
                    lua_State *L,
                    Qx(QLA_D,_DiracPropagator) *r,
                    int mu, int sign,
                    Qx(QLA_D,_DiracFermion) **a)
#endif
{
#if QNc == 'N'
    typedef QLA_DN_DiracFermion(nc, Vtype);
    typedef QLA_DN_HalfFermion(nc, Htype);
#else
    typedef Qx(QLA_D,_DiracFermion) Vtype;
    typedef Qx(QLA_D,_HalfFermion) Htype;
#endif
    int count = QDP_sites_on_node;
    int k, ic, is;

    if (mu == 0 || mu == 2) sign = 1 - sign; /* AAA gamma basis dependent */

    for (k = 0; k < count; k++) {
        for (ic = 0; ic < nc; ic++) {
            for (is = 0; is < QDP_Ns; is++) {
                int jc, js;
                Vtype fk;
                Htype hk;
                for (jc = 0; jc < nc; jc++) {
                    for (js = 0; js < QDP_Ns/2; js++) {
                        QLA_c_eq_c(QLA_elem_H(hk, jc, js),
                                   QLA_elem_D(*a[PIDX(jc,js)], ic, is));
                    }
                }
#if QNc == 'N'
                Qx(QLA_D,_D_eq_sprecon_H)(nc, &fk, &hk, mu, sign);
#else
                Qx(QLA_D,_D_eq_sprecon_H)(&fk, &hk, mu, sign);
#endif
                for (jc = 0; jc < nc; jc++) {
                    for (js = 0; js < QDP_Ns; js++) {
                        QLA_c_eq_c(QLA_elem_P(*r, ic, is, jc, js),
                                   QLA_elem_D(fk, jc, js));
                    }
                }
            }
        }
        r++;
        for (ic = 0; ic < nc; ic++)
            for (is = 0; is < QDP_Ns/2; is++)
                a[PIDX(ic,is)]++;
    }    
}

static int
Qs(q_P_reconstruct)(lua_State *L, Qs(X_reconstruct) *op, int nc)
{
    const char *sign = luaL_checkstring(L, 1);
    int mu = qlua_checkgammaindex(L, 2);
    Qs(mLatDirProp) *r;
    int isign = 0;
    int ic, is;
    mLattice *S = NULL;
    int Sidx = 0;
    Qs(mLatDirFerm) *a[nc * QDP_Ns/2];
#if QNc == 'N'
    QLA_DN_DiracPropagator(nc, (*rr));
    QLA_DN_DiracFermion(nc, (*aa[nc * QDP_Ns/2]));
#else
    Qx(QLA_D,_DiracPropagator) *rr;
    Qx(QLA_D,_DiracFermion) *aa[nc * QDP_Ns/2];
#endif
    
    luaL_checktype(L, 3, LUA_TTABLE);
    for (ic = 0; ic < nc; ic++) {
        lua_pushnumber(L, ic + 1);
        lua_gettable(L, 3);
        int tic_idx = lua_gettop(L);
        qlua_checktable(L, -1, "projected Prop" Qcolors " (%d,*)", ic);
        for (is = 0; is < QDP_Ns/2; is++) {
            lua_pushnumber(L, is + 1);
            lua_gettable(L, tic_idx);
            a[PIDX(ic, is)] = Qs(qlua_checkLatDirFerm)(L, -1, S, nc);
            if (S == 0) {
                S = qlua_ObjLattice(L, -1);
                Sidx = lua_gettop(L);
            }
        }
    }

    if (strcmp(sign, "plus") == 0)
        isign = 1;
    else if (strcmp(sign, "minus") == 0)
        isign = 0;
    else
        return luaL_error(L, "bad sign parameter");

    r = Qs(qlua_newLatDirProp)(L, Sidx, nc);

    if (mu == 5) mu = 4; /* [sic] -- funny numbering in QLA's spproj/sprecon */

    CALL_QDP(L);
    /* ZZZ will break if any of a[][] is aliased */
    for (ic = 0; ic < nc; ic++) {
        for (is = 0; is < QDP_Ns/2; is++) {
            int idx_cs = PIDX(ic, is);
            aa[idx_cs] = Qx(QDP_D,_expose_D)(a[idx_cs]->ptr);
        }
    }
    rr = Qx(QDP_D,_expose_P)(r->ptr);

    op(nc, L, rr, mu, isign, aa);

    Qx(QDP_D,_reset_P)(r->ptr);
    for (ic = 0; ic < nc; ic++) {
        for (is = 0; is < QDP_Ns/2; is++) {
            Qx(QDP_D,_reset_D)(a[PIDX(ic, is)]->ptr);
        }
    }
    return 1;
}

#undef PIDX
#undef QNc
#undef Qcolors
#undef Qs
#undef Qx
#undef QC
