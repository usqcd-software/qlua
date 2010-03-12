#ifndef MARK_C55CC38A_5878_4666_B0BE_088EEFB573A3
#define MARK_C55CC38A_5878_4666_B0BE_088EEFB573A3

typedef struct {
    QDP_DiracFermion *ptr;
} mLatDirFerm;

extern const char mtnLatDirFerm[];

int init_latdirferm(lua_State *L);
int fini_latdirferm(lua_State *L);

mLatDirFerm *qlua_checkLatDirFerm(lua_State *L, int idx);
mLatDirFerm *qlua_newLatDirFerm(lua_State *L);

int q_D_gaussian(lua_State *L);

#endif /* !defined(MARK_C55CC38A_5878_4666_B0BE_088EEFB573A3) */
