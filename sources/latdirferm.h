#ifndef MARK_C55CC38A_5878_4666_B0BE_088EEFB573A3
#define MARK_C55CC38A_5878_4666_B0BE_088EEFB573A3

#if USE_Nc2
typedef struct {
    QDP_D2_DiracFermion *ptr;
} mLatDirFerm2;

mLatDirFerm2 *qlua_checkLatDirFerm2(lua_State *L, int idx, mLattice *S, int nc);
mLatDirFerm2 *qlua_newLatDirFerm2(lua_State *L, int Sidx, int nc);
mLatDirFerm2 *qlua_newZeroLatDirFerm2(lua_State *L, int Sidx, int nc);
#endif

#if USE_Nc3
typedef struct {
    QDP_D3_DiracFermion *ptr;
} mLatDirFerm3;

mLatDirFerm3 *qlua_checkLatDirFerm3(lua_State *L, int idx, mLattice *S, int nc);
mLatDirFerm3 *qlua_newLatDirFerm3(lua_State *L, int Sidx, int nc);
mLatDirFerm3 *qlua_newZeroLatDirFerm3(lua_State *L, int Sidx, int nc);
#endif

#if USE_NcN
typedef struct {
    int nc;
    QDP_DN_DiracFermion *ptr;
} mLatDirFermN;

mLatDirFermN *qlua_checkLatDirFermN(lua_State *L, int idx, mLattice *S, int nc);
mLatDirFermN *qlua_newLatDirFermN(lua_State *L, int Sidx, int nc);
mLatDirFermN *qlua_newZeroLatDirFermN(lua_State *L, int Sidx, int nc);
#endif

int init_latdirferm(lua_State *L);
void fini_latdirferm(void);

int q_D_gaussian(lua_State *L);
int q_D_gaussian_N(lua_State *L);

#endif /* !defined(MARK_C55CC38A_5878_4666_B0BE_088EEFB573A3) */
