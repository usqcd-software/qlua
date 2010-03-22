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
#if QNc == 'N'
    Qx(QDP_D,_D_eq_zero)(QC(f), r->ptr, *S->qss);
#else
    Qx(QDP_D,_D_eq_zero)(r->ptr, *S->qss);
#endif
    for (i = 0; i < 16; i++) {
        switch (m->g[i].t) {
        case qG_z: continue;
        case qG_p:
#if QNc == 'N'
            Qx(QDP_D,_D_eq_gamma_times_D)(QC(f), mf->ptr, f->ptr, i, *S->qss);
            Qx(QDP_D,_D_peq_D)(QC(f), r->ptr, mf->ptr, *S->qss);
#else
            Qx(QDP_D,_D_eq_gamma_times_D)(mf->ptr, f->ptr, i, *S->qss);
            Qx(QDP_D,_D_peq_D)(r->ptr, mf->ptr, *S->qss);
#endif
            break;
        case qG_m:
#if QNc == 'N'
            Qx(QDP_D,_D_eq_gamma_times_D)(QC(f), mf->ptr, f->ptr, i, *S->qss);
            Qx(QDP_D,_D_meq_D)(QC(f), r->ptr, mf->ptr, *S->qss);
#else
            Qx(QDP_D,_D_eq_gamma_times_D)(mf->ptr, f->ptr, i, *S->qss);
            Qx(QDP_D,_D_meq_D)(r->ptr, mf->ptr, *S->qss);
#endif
            break;
        case qG_r:
#if QNc == 'N'
            Qx(QDP_D,_D_eq_gamma_times_D)(QC(f), mf->ptr, f->ptr, i, *S->qss);
            Qx(QDP_D,_D_peq_r_times_D)(QC(f), r->ptr, &m->g[i].r, mf->ptr,
                                       *S->qss);
#else
            Qx(QDP_D,_D_eq_gamma_times_D)(mf->ptr, f->ptr, i, *S->qss);
            Qx(QDP_D,_D_peq_r_times_D)(r->ptr, &m->g[i].r, mf->ptr, *S->qss);
#endif
            break;
        case qG_c:
#if QNc == 'N'
            Qx(QDP_D,_D_eq_gamma_times_D)(QC(f), mf->ptr, f->ptr, i, *S->qss);
            Qx(QDP_D,_D_peq_c_times_D)(QC(f), r->ptr, &m->g[i].c, mf->ptr,
                                       *S->qss);
#else
            Qx(QDP_D,_D_eq_gamma_times_D)(mf->ptr, f->ptr, i, *S->qss);
            Qx(QDP_D,_D_peq_c_times_D)(r->ptr, &m->g[i].c, mf->ptr, *S->qss);
#endif
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
#if QNc == 'N'
    Qx(QDP_D,_P_eq_zero)(QC(f), r->ptr, *S->qss);
#else
    Qx(QDP_D,_P_eq_zero)(r->ptr, *S->qss);
#endif
    for (i = 0; i < 16; i++) {
        switch (m->g[i].t) {
        case qG_z: continue;
        case qG_p:
#if QNc == 'N'
            Qx(QDP_D,_P_eq_gamma_times_P)(QC(f), mf->ptr, f->ptr, i,
                                          *S->qss);
            Qx(QDP_D,_P_peq_P)(QC(f), r->ptr, mf->ptr, *S->qss);
#else
            Qx(QDP_D,_P_eq_gamma_times_P)(mf->ptr, f->ptr, i,
                                          *S->qss);
            Qx(QDP_D,_P_peq_P)(r->ptr, mf->ptr, *S->qss);
#endif
            break;
        case qG_m:
#if QNc == 'N'
            Qx(QDP_D,_P_eq_gamma_times_P)(QC(f), mf->ptr, f->ptr, i,
                                          *S->qss);
            Qx(QDP_D,_P_meq_P)(QC(f), r->ptr, mf->ptr, *S->qss);
#else
            Qx(QDP_D,_P_eq_gamma_times_P)(mf->ptr, f->ptr, i,
                                          *S->qss);
            Qx(QDP_D,_P_meq_P)(r->ptr, mf->ptr, *S->qss);
#endif
            break;
        case qG_r:
#if QNc == 'N'
            Qx(QDP_D,_P_eq_gamma_times_P)(QC(f), mf->ptr, f->ptr, i,
                                          *S->qss);
            Qx(QDP_D,_P_peq_r_times_P)(QC(f), r->ptr, &m->g[i].r, mf->ptr,
                                       *S->qss);
#else
            Qx(QDP_D,_P_eq_gamma_times_P)(mf->ptr, f->ptr, i,
                                          *S->qss);
            Qx(QDP_D,_P_peq_r_times_P)(r->ptr, &m->g[i].r, mf->ptr,
                                       *S->qss);
#endif
            break;
        case qG_c:
#if QNc == 'N'
               Qx(QDP_D,_P_eq_gamma_times_P)(QC(f), mf->ptr, f->ptr, i,
                                             *S->qss);
               Qx(QDP_D,_P_peq_c_times_P)(QC(f), r->ptr, &m->g[i].c, mf->ptr,
                                          *S->qss);
#else
               Qx(QDP_D,_P_eq_gamma_times_P)(mf->ptr, f->ptr, i,
                                             *S->qss);
               Qx(QDP_D,_P_peq_c_times_P)(r->ptr, &m->g[i].c, mf->ptr,
                                          *S->qss);
#endif
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
#if QNc == 'N'
    Qx(QDP_D,_P_eq_zero)(QC(f), r->ptr, *S->qss);
#else
    Qx(QDP_D,_P_eq_zero)(r->ptr, *S->qss);
#endif
    for (i = 0; i < 16; i++) {
        switch (m->g[i].t) {
        case qG_z: continue;
        case qG_p:
#if QNc == 'N'
            Qx(QDP_D,_P_eq_P_times_gamma)(QC(f), mf->ptr, f->ptr, i,
                                          *S->qss);
            Qx(QDP_D,_P_peq_P)(QC(f), r->ptr, mf->ptr,
                               *S->qss);
#else
            Qx(QDP_D,_P_eq_P_times_gamma)(mf->ptr, f->ptr, i, *S->qss);
            Qx(QDP_D,_P_peq_P)(r->ptr, mf->ptr, *S->qss);
#endif
            break;
        case qG_m:
#if QNc == 'N'
            Qx(QDP_D,_P_eq_P_times_gamma)(QC(f), mf->ptr, f->ptr, i,
                                          *S->qss);
            Qx(QDP_D,_P_meq_P)(QC(f), r->ptr, mf->ptr,
                               *S->qss);
#else
            Qx(QDP_D,_P_eq_P_times_gamma)(mf->ptr, f->ptr, i, *S->qss);
            Qx(QDP_D,_P_meq_P)(r->ptr, mf->ptr, *S->qss);
#endif
            break;
        case qG_r:
#if QNc == 'N'
            Qx(QDP_D,_P_eq_P_times_gamma)(QC(f), mf->ptr, f->ptr, i,
                                          *S->qss);
            Qx(QDP_D,_P_peq_r_times_P)(QC(f), r->ptr, &m->g[i].r, mf->ptr,
                                       *S->qss);
#else
            Qx(QDP_D,_P_eq_P_times_gamma)(mf->ptr, f->ptr, i, *S->qss);
            Qx(QDP_D,_P_peq_r_times_P)(r->ptr, &m->g[i].r, mf->ptr, *S->qss);
#endif
            break;
        case qG_c:
#if QNc == 'N'
            Qx(QDP_D,_P_eq_P_times_gamma)(QC(f), mf->ptr, f->ptr, i,
                                          *S->qss);
            Qx(QDP_D,_P_peq_c_times_P)(QC(f), r->ptr, &m->g[i].c, mf->ptr,
                                       *S->qss);
#else
            Qx(QDP_D,_P_eq_P_times_gamma)(mf->ptr, f->ptr, i, *S->qss);
            Qx(QDP_D,_P_peq_c_times_P)(r->ptr, &m->g[i].c, mf->ptr, *S->qss);
#endif
            break;
        }
    }
    return 1;
}

typedef void Qs(X_project)(lua_State *L,
                           Qx(QLA_D,_DiracFermion) **r,
                           int mu, int sign,
                           Qx(QLA_D,_DiracPropagator) *f,
                           int nc);

static void
Qs(q_P_left_proj)(lua_State *L,
                  Qx(QLA_D,_DiracFermion) **r,
                  int mu,
                  int sign,
                  Qx(QLA_D,_DiracPropagator) *f,
                  int nc)
{
    int count = QDP_sites_on_node;
    int k, ic, is;

    for (k = 0; k < count; k++) {
        for (ic = 0; ic < nc; ic++) {
            for (is = 0; is < QDP_Ns; is++) {
                int jc, js;
                Qx(QLA_D,_DiracFermion) fk;
                Qx(QLA_D,_HalfFermion) hk;
                for (jc = 0; jc < nc; jc++) {
                    for (js = 0; js < QDP_Ns; js++) {
                        QLA_c_eq_c(Qx(QLA_D,_elem_D)(fk, jc, js),
                                   Qx(QLA_D,_elem_P)(*f, jc, js, ic, is));
                    }
                }
#if QNc == 'N'
                Qx(QLA_D,_H_eq_spproj_D)(nc, &hk, &fk, mu, sign);
#else
                Qx(QLA_D,_H_eq_spproj_D)(&hk, &fk, mu, sign);
#endif
                for (jc = 0; jc < nc; jc++) {
                    for (js = 0; js < QDP_Ns/2; js++) {
                        QLA_c_eq_c(Qx(QLA_D,_elem_D)(*r[PIDX(jc, js)], ic, is),
                                   Qx(QLA_D,_elem_H)(hk, jc, js));
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
Qs(q_P_right_proj)(lua_State *L,
                   Qx(QLA_D,_DiracFermion) **r,
                   int mu,
                   int sign,
                   Qx(QLA_D,_DiracPropagator) *f,
                   int nc)
{
    int count = QDP_sites_on_node;
    int k, ic, is;

    if (mu == 0 || mu == 2) sign = 1 - sign; /* AAA gamma basis dependent */

    for (k = 0; k < count; k++) {
        for (ic = 0; ic < nc; ic++) {
            for (is = 0; is < QDP_Ns; is++) {
                int jc, js;
                Qx(QLA_D,_DiracFermion) fk;
                Qx(QLA_D,_HalfFermion) hk;
                for (jc = 0; jc < nc; jc++) {
                    for (js = 0; js < QDP_Ns; js++) {
                        QLA_c_eq_c(Qx(QLA_D,_elem_D)(fk, jc, js),
                                   Qx(QLA_D,_elem_P)(*f, ic, is, jc, js));
                    }
                }
#if QNc == 'N'
                Qx(QLA_D,_H_eq_spproj_D)(nc, &hk, &fk, mu, sign);
#else
                Qx(QLA_D,_H_eq_spproj_D)(&hk, &fk, mu, sign);
#endif
                for (jc = 0; jc < nc; jc++) {
                    for (js = 0; js < QDP_Ns/2; js++) {
                        QLA_c_eq_c(Qx(QLA_D,_elem_D)(*r[PIDX(jc, js)], ic, is),
                                   Qx(QLA_D,_elem_H)(hk, jc, js));
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
    Qx(QLA_D,_DiracPropagator) *ff;
    Qs(mLatDirFerm) **r;
    Qx(QLA_D,_DiracFermion) **rr;

    qlua_ObjLattice(L, 3);
    Sidx = lua_gettop(L);
    r = qlua_malloc(L, QC(f) * QDP_Ns/2 * sizeof (Qs(mLatDirFerm) *));
    rr = qlua_malloc(L, QC(f) * QDP_Ns/2 * sizeof (Qx(QLA_D,_DiracFermion) *));

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
#if QNc == 'N'
            rr[idx_cs] = Qx(QDP_D,_expose_D)(QC(f), r[idx_cs]->ptr);
#else
            rr[idx_cs] = Qx(QDP_D,_expose_D)(r[idx_cs]->ptr);
#endif
        }
    }
#if QNc == 'N'
    ff = Qx(QDP_D,_expose_P)(QC(f), f->ptr);
#else
    ff = Qx(QDP_D,_expose_P)(f->ptr);
#endif

    op(L, rr, mu, isign, ff, QC(f));

#if QNc == 'N'
    Qx(QDP_D,_reset_P)(QC(f), f->ptr);
#else
    Qx(QDP_D,_reset_P)(f->ptr);
#endif
    for (ic = 0; ic < QC(f); ic++) {
        for (is = 0; is < QDP_Ns/2; is++)  {
            int idx_cs = PIDX(ic, is);
#if QNc == 'N'
            Qx(QDP_D,_reset_D)(QC(f), r[idx_cs]->ptr);
#else
            Qx(QDP_D,_reset_D)(r[idx_cs]->ptr);
#endif
        }
    }
    qlua_free(L, r);
    qlua_free(L, rr);

    return 1;
}

typedef void Qs(X_reconstruct)(lua_State *L,
                               Qx(QLA_D,_DiracPropagator) *r,
                               int mu, int sign,
                               Qx(QLA_D,_DiracFermion) **a,
                               int nc);

static void
Qs(q_P_left_recon)(lua_State *L,
                   Qx(QLA_D,_DiracPropagator) *r,
                   int mu, int sign,
                   Qx(QLA_D,_DiracFermion) **a,
                   int nc)
{
    int count = QDP_sites_on_node;
    int k, ic, is;

    for (k = 0; k < count; k++) {
        for (ic = 0; ic < nc; ic++) {
            for (is = 0; is < QDP_Ns; is++) {
                int jc, js;
                Qx(QLA_D,_DiracFermion) fk;
                Qx(QLA_D,_HalfFermion) hk;
                for (jc = 0; jc < nc; jc++) {
                    for (js = 0; js < QDP_Ns/2; js++) {
                        QLA_c_eq_c(Qx(QLA_D,_elem_H)(hk, jc, js),
                                   Qx(QLA_D,_elem_D)(*a[PIDX(jc,js)], ic, is));
                    }
                }
#if QNc == 'N'
                Qx(QLA_D,_D_eq_sprecon_H)(nc, &fk, &hk, mu, sign);
#else
                Qx(QLA_D,_D_eq_sprecon_H)(&fk, &hk, mu, sign);
#endif
                for (jc = 0; jc < QDP_Nc; jc++) {
                    for (js = 0; js < QDP_Ns; js++) {
                        QLA_c_eq_c(Qx(QLA_D,_elem_P)(*r, jc, js, ic, is),
                                   Qx(QLA_D,_elem_D)(fk, jc, js));
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

static void
Qs(q_P_right_recon)(lua_State *L,
                    Qx(QLA_D,_DiracPropagator) *r,
                    int mu, int sign,
                    Qx(QLA_D,_DiracFermion) **a,
                    int nc)
{
    int count = QDP_sites_on_node;
    int k, ic, is;

    if (mu == 0 || mu == 2) sign = 1 - sign; /* AAA gamma basis dependent */

    for (k = 0; k < count; k++) {
        for (ic = 0; ic < nc; ic++) {
            for (is = 0; is < QDP_Ns; is++) {
                int jc, js;
                Qx(QLA_D,_DiracFermion) fk;
                Qx(QLA_D,_HalfFermion) hk;
                for (jc = 0; jc < QDP_Nc; jc++) {
                    for (js = 0; js < QDP_Ns/2; js++) {
                        QLA_c_eq_c(Qx(QLA_D,_elem_H)(hk, jc, js),
                                   Qx(QLA_D,_elem_D)(*a[PIDX(jc,js)], ic, is));
                    }
                }
#if QNc == 'N'
                Qx(QLA_D,_D_eq_sprecon_H)(nc, &fk, &hk, mu, sign);
#else
                Qx(QLA_D,_D_eq_sprecon_H)(&fk, &hk, mu, sign);
#endif
                for (jc = 0; jc < nc; jc++) {
                    for (js = 0; js < QDP_Ns; js++) {
                        QLA_c_eq_c(Qx(QLA_D,_elem_P)(*r, ic, is, jc, js),
                                   Qx(QLA_D,_elem_D)(fk, jc, js));
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
    Qs(mLatDirFerm) **a;
    Qx(QLA_D,_DiracFermion) **aa;
    Qx(QLA_D,_DiracPropagator) *rr;
    mLattice *S = NULL;
    int Sidx;
    
    a = qlua_malloc(L, nc * QDP_Ns/2 * sizeof (Qs(mLatDirFerm) *));
    aa = qlua_malloc(L, nc * QDP_Ns/2 * sizeof (Qx(QLA_D,_DiracFermion) *));
    luaL_checktype(L, 3, LUA_TTABLE);
    for (ic = 0; ic < nc; ic++) {
        lua_pushnumber(L, ic + 1);
        lua_gettable(L, 3);
        qlua_checktable(L, -1, "projected Prop" Qcolors " (%d,*)", ic);
        for (is = 0; is < QDP_Ns/2; is++) {
            lua_pushnumber(L, is + 1);
            lua_gettable(L, -2);
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
#if QNc == 'N'
            aa[idx_cs] = Qx(QDP_D,_expose_D)(nc, a[idx_cs]->ptr);
#else
            aa[idx_cs] = Qx(QDP_D,_expose_D)(a[idx_cs]->ptr);
#endif
        }
    }
#if QNc == 'N'
    rr = Qx(QDP_D,_expose_P)(nc, r->ptr);
#else
    rr = Qx(QDP_D,_expose_P)(r->ptr);
#endif

    op(L, rr, mu, isign, aa, nc);

#if QNc == 'N'
    Qx(QDP_D,_reset_P)(nc, r->ptr);
#else
    Qx(QDP_D,_reset_P)(r->ptr);
#endif
    for (ic = 0; ic < nc; ic++) {
        for (is = 0; is < QDP_Ns/2; is++) {
#if QNc == 'N'
            Qx(QDP_D,_reset_D)(nc, a[PIDX(ic, is)]->ptr);
#else
            Qx(QDP_D,_reset_D)(a[PIDX(ic, is)]->ptr);
#endif
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
