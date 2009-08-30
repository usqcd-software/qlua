#ifndef MARK_1487E8E3_E4CD_4CEB_80B2_5BA83A4156D4
#define MARK_1487E8E3_E4CD_4CEB_80B2_5BA83A4156D4

typedef struct {
    QDP_Int *ptr;
} mLatInt;

extern const char mtnLatInt[];

int init_latint(lua_State *L);
int fini_latint(lua_State *L);

mLatInt *qlua_checkLatInt(lua_State *L, int idx);
mLatInt *qlua_newLatInt(lua_State *L);

int qlua_index(lua_State *L, int idx, const char *name, int mv); /* k or -1 */
int qlua_checkindex(lua_State *L, int idx, const char *name, int mv);
int qlua_diracindex(lua_State *L, int idx, int mv); /* k or -1 */
int qlua_checkdiracindex(lua_State *L, int idx, int mv);
int qlua_colorindex(lua_State *L, int idx, int mv); /* k or -1 */
int qlua_checkcolorindex(lua_State *L, int idx, int mv);
int *qlua_latcoord(lua_State *L, int idx); /* lc[] or NULL */
int *qlua_checklatcoord(lua_State *L, int idx);
int qlua_checkleftindex(lua_State *L, int idx, int mv);
int qlua_leftindex(lua_State *L, int idx, int mv);
int qlua_checkrightindex(lua_State *L, int idx, int mv);
int qlua_rightindex(lua_State *L, int idx, int mv);
QDP_Shift qlua_checkShift(lua_State *L, int idx);
QDP_ShiftDir qlua_checkShiftDir(lua_State *L, int idx);

int q_I_add_I(lua_State *L);
int q_I_sub_I(lua_State *L);
int q_I_mul_I(lua_State *L);
int q_i_mul_I(lua_State *L);
int q_I_mul_i(lua_State *L);
int q_I_div_I(lua_State *L);

#endif /* !defined(MARK_1487E8E3_E4CD_4CEB_80B2_5BA83A4156D4) */
