#ifndef MARK_1EF46601_F4FC_44FF_8701_E58341017024
#define MARK_1EF46601_F4FC_44FF_8701_E58341017024

typedef struct {
    QDP_D_Complex *ptr;
} mLatComplex;

extern const char *mtnLatComplex;

int init_latcomplex(lua_State *L);
int fini_latcomplex(lua_State *L);

mLatComplex *qlua_checkLatComplex(lua_State *L, int idx);
mLatComplex *qlua_newLatComplex(lua_State *L);

int q_C_gaussian(lua_State *L);

int q_C_dot(lua_State *L);
int q_C_add_C(lua_State *L);
int q_C_sub_C(lua_State *L);
int q_C_mul_C(lua_State *L);
int q_R_mul_C(lua_State *L);
int q_C_mul_R(lua_State *L);
int q_c_mul_C(lua_State *L);
int q_C_mul_c(lua_State *L);
int q_r_mul_C(lua_State *L);
int q_C_mul_r(lua_State *L);
int q_C_div_C(lua_State *L);


#endif /* !defined(MARK_1EF46601_F4FC_44FF_8701_E58341017024) */
