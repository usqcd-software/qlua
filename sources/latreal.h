#ifndef MARK_DA4C72A1_DB3E_44A4_931E_369DDFF0D309
#define MARK_DA4C72A1_DB3E_44A4_931E_369DDFF0D309

typedef struct {
    QDP_Real *ptr;
} mLatReal;

extern const char *mtnLatReal;

int init_latreal(lua_State *L);
int fini_latreal(lua_State *L);

mLatReal *qlua_checkLatReal(lua_State *L, int idx);
mLatReal *qlua_newLatReal(lua_State *L);

int q_R_random(lua_State *L);
int q_R_gaussian(lua_State *L);

int q_R_add_R(lua_State *L);
int q_R_sub_R(lua_State *L);
int q_R_mul_R(lua_State *L);
int q_r_mul_R(lua_State *L);
int q_R_mul_r(lua_State *L);
int q_R_div_R(lua_State *L);

#endif /* !defined(MARK_DA4C72A1_DB3E_44A4_931E_369DDFF0D309) */
