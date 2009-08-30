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

int q_D_dot(lua_State *L);
int q_D_add_D(lua_State *L);
int q_D_sub_D(lua_State *L);
int q_r_mul_D(lua_State *L);
int q_D_mul_r(lua_State *L);
int q_c_mul_D(lua_State *L);
int q_D_mul_c(lua_State *L);
int q_M_mul_D(lua_State *L);

#endif /* !defined(MARK_C55CC38A_5878_4666_B0BE_088EEFB573A3) */
