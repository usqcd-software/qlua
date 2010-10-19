#ifndef MARK_1487E8E3_E4CD_4CEB_80B2_5BA83A4156D4
#define MARK_1487E8E3_E4CD_4CEB_80B2_5BA83A4156D4

typedef struct {
    QDP_Int *ptr;
} mLatInt;

int init_latint(lua_State *L);
int fini_latint(lua_State *L);

mLatInt *qlua_checkLatInt(lua_State *L, int idx, mLattice *S);
mLatInt *qlua_newLatInt(lua_State *L, int S_idx);

#endif /* !defined(MARK_1487E8E3_E4CD_4CEB_80B2_5BA83A4156D4) */
