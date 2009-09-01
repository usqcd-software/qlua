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

QDP_Shift qlua_checkShift(lua_State *L, int idx);
QDP_ShiftDir qlua_checkShiftDir(lua_State *L, int idx);

#endif /* !defined(MARK_1487E8E3_E4CD_4CEB_80B2_5BA83A4156D4) */
