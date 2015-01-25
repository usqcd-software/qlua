#ifndef MARK_1487E8E3_E4CD_4CEB_80B2_5BA83A4156D4
#define MARK_1487E8E3_E4CD_4CEB_80B2_5BA83A4156D4

typedef struct {
    QDP_Int *ptr;
} mLatInt;

int init_latint(lua_State *L);
void fini_latint(void);

mLatInt *qlua_checkLatInt(lua_State *L, int idx, mLattice *S);
mLatInt *qlua_newLatInt(lua_State *L, int S_idx);
mLatInt *qlua_newZeroLatInt(lua_State *L, int S_idx);

/* NB: be careful not to drop the result object from memory */
mLatInt *qlua_tabkey_LatInt(lua_State *L, int idx, const char *key, mLattice *S);
mLatInt *qlua_tabidx_LatInt(lua_State *L, int idx, int subidx, mLattice *S);

#endif /* !defined(MARK_1487E8E3_E4CD_4CEB_80B2_5BA83A4156D4) */
