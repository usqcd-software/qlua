#ifndef MARK_DA4C72A1_DB3E_44A4_931E_369DDFF0D309
#define MARK_DA4C72A1_DB3E_44A4_931E_369DDFF0D309

typedef struct {
    QDP_Real *ptr;
} mLatReal;

int init_latreal(lua_State *L);
int fini_latreal(lua_State *L);

mLatReal *qlua_checkLatReal(lua_State *L, int idx, mLattice *S);
mLatReal *qlua_newLatReal(lua_State *L, int Sidx);
mLatReal *qlua_newZeroLatReal(lua_State *L, int Sidx);

int q_R_random(lua_State *L);
int q_R_gaussian(lua_State *L);

#endif /* !defined(MARK_DA4C72A1_DB3E_44A4_931E_369DDFF0D309) */
