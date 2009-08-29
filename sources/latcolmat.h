#ifndef MARK_E53A3544_2065_4CA4_888F_72B524D3269B
#define MARK_E53A3544_2065_4CA4_888F_72B524D3269B

typedef struct {
    QDP_ColorMatrix *ptr;
} mLatColMat;

extern const char *mtnLatColMat;

int init_latcolmat(lua_State *L);
int fini_latcolmat(lua_State *L);

mLatColMat *qlua_checkLatColMat(lua_State *L, int idx);
mLatColMat *qlua_newLatColMat(lua_State *L);

int q_M_gaussian(lua_State *L);

int q_M_dot(lua_State *L);
int q_M_add_M(lua_State *L);
int q_M_sub_M(lua_State *L);
int q_r_mul_M(lua_State *L);
int q_M_mul_r(lua_State *L);
int q_c_mul_M(lua_State *L);
int q_M_mul_c(lua_State *L);
int q_M_mul_V(lua_State *L);
int q_M_mul_M(lua_State *L);


#endif /* !defined(MARK_E53A3544_2065_4CA4_888F_72B524D3269B) */
