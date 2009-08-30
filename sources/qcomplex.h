#ifndef MARK_91468699_C123_4148_86B5_3DAB4308E77E
#define MARK_91468699_C123_4148_86B5_3DAB4308E77E
#include <qla.h>

extern const char mtnComplex[];

int init_complex(lua_State *L);
int fini_complex(lua_State *L);

QLA_Complex *qlua_checkComplex(lua_State *L, int idx);
QLA_Complex *qlua_newComplex(lua_State *L);

/* additions */
int q_r_add_c(lua_State *L);
int q_c_add_r(lua_State *L);
int q_c_add_c(lua_State *L);

/* subtractions */
int q_r_sub_c(lua_State *L);
int q_c_sub_r(lua_State *L);
int q_c_sub_c(lua_State *L);

/* multiplications */
int q_r_mul_c(lua_State *L);
int q_c_mul_r(lua_State *L);
int q_c_mul_c(lua_State *L);

/* divisions */
int q_r_div_c(lua_State *L);
int q_c_div_r(lua_State *L);
int q_c_div_c(lua_State *L);


#endif /* !defined(MARK_91468699_C123_4148_86B5_3DAB4308E77E) */
