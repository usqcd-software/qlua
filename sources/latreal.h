#ifndef MARK_DA4C72A1_DB3E_44A4_931E_369DDFF0D309
#define MARK_DA4C72A1_DB3E_44A4_931E_369DDFF0D309

typedef struct {
    QDP_D_Real *ptr;
} mLatRealD;

typedef struct {
    QDP_F_Real *ptr;
} mLatRealF;

int init_latreal(lua_State *L);
int fini_latreal(lua_State *L);

mLatRealD *qlua_checkLatRealD(lua_State *L, int idx, mLattice *S);
mLatRealF *qlua_checkLatRealF(lua_State *L, int idx, mLattice *S);
mLatRealD *qlua_newLatRealD(lua_State *L, int Sidx);
mLatRealF *qlua_newLatRealF(lua_State *L, int Sidx);

int q_R_random(lua_State *L);
int q_RD_random(lua_State *L);
int q_RF_random(lua_State *L);
int q_R_gaussian(lua_State *L);
int q_RD_gaussian(lua_State *L);
int q_RF_gaussian(lua_State *L);

#endif /* !defined(MARK_DA4C72A1_DB3E_44A4_931E_369DDFF0D309) */
