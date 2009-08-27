#ifndef MARK_1487E8E3_E4CD_4CEB_80B2_5BA83A4156D4
#define MARK_1487E8E3_E4CD_4CEB_80B2_5BA83A4156D4

typedef struct {
    QDP_Int *ptr;
} mLatInt;

extern const char *mtnLatInt;

int init_latint(lua_State *L);
int fini_latint(lua_State *L);

mLatInt *q_checkLatInt(lua_State *L, int idx);
mLatInt *q_newLatInt(lua_State *L);

/* copies */
int q_I_eq_I(lua_State *L);

/* additions */
int q_I_add_I(lua_State *L);

/* subtractions */
int q_I_sub_I(lua_State *L);

/* multiplications */
int q_I_mul_I(lua_State *L);
int q_i_mul_I(lua_State *L);
int q_I_mul_i(lua_State *L);

/* divisions */
int q_I_div_I(lua_State *L);

#endif /* !defined(MARK_1487E8E3_E4CD_4CEB_80B2_5BA83A4156D4) */
