#ifndef MARK_5F7EAA43_B7D7_4A36_B574_D9707797522B
#define MARK_5F7EAA43_B7D7_4A36_B574_D9707797522B

#if USE_Nc2
typedef struct {
    QDP_D2_DiracPropagator *ptr;
} mLatDirProp2;

mLatDirProp2 *qlua_checkLatDirProp2(lua_State *L, int idx, mLattice *S, int nc);
mLatDirProp2 *qlua_newLatDirProp2(lua_State *L, int Sidx, int nc);
mLatDirProp2 *qlua_newZeroLatDirProp2(lua_State *L, int Sidx, int nc);
#endif

#if USE_Nc3
typedef struct {
    QDP_D3_DiracPropagator *ptr;
} mLatDirProp3;

mLatDirProp3 *qlua_checkLatDirProp3(lua_State *L, int idx, mLattice *S, int nc);
mLatDirProp3 *qlua_newLatDirProp3(lua_State *L, int Sidx, int nc);
mLatDirProp3 *qlua_newZeroLatDirProp3(lua_State *L, int Sidx, int nc);
#endif

#if USE_NcN
typedef struct {
    int nc;
    QDP_DN_DiracPropagator *ptr;
} mLatDirPropN;

mLatDirPropN *qlua_checkLatDirPropN(lua_State *L, int idx, mLattice *S, int nc);
mLatDirPropN *qlua_newLatDirPropN(lua_State *L, int Sidx, int nc);
mLatDirPropN *qlua_newZeroLatDirPropN(lua_State *L, int Sidx, int nc);
#endif

int init_latdirprop(lua_State *L);
int fini_latdirprop(lua_State *L);

int q_P_gaussian(lua_State *L);
int q_P_gaussian_N(lua_State *L);

#endif /* !defined(MARK_5F7EAA43_B7D7_4A36_B574_D9707797522B) */
