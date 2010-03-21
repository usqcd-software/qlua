#ifndef MARK_C55CC38A_5878_4666_B0BE_088EEFB573A3
#define MARK_C55CC38A_5878_4666_B0BE_088EEFB573A3

typedef struct {
    QDP_D2_DiracFermion *ptr;
} mLatDirFerm2;

typedef struct {
    QDP_D3_DiracFermion *ptr;
} mLatDirFerm3;

typedef struct {
    int nc;
    QDP_DN_DiracFermion *ptr;
} mLatDirFermN;

int init_latdirferm(lua_State *L);
int fini_latdirferm(lua_State *L);

mLatDirFerm2 *qlua_checkLatDirFerm2(lua_State *L, int idx, mLattice *S, int nc);
mLatDirFerm3 *qlua_checkLatDirFerm3(lua_State *L, int idx, mLattice *S, int nc);
mLatDirFermN *qlua_checkLatDirFermN(lua_State *L, int idx, mLattice *S, int nc);
mLatDirFerm2 *qlua_newLatDirFerm2(lua_State *L, int Sidx, int nc);
mLatDirFerm3 *qlua_newLatDirFerm3(lua_State *L, int Sidx, int nc);
mLatDirFermN *qlua_newLatDirFermN(lua_State *L, int Sidx, int nc);

int q_D_gaussian(lua_State *L);
int q_D_gaussian_N(lua_State *L);

#endif /* !defined(MARK_C55CC38A_5878_4666_B0BE_088EEFB573A3) */
