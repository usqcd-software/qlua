#ifndef MARK_C8FCBCFB_8280_4F1D_9159_8E1C3F95E045
#define MARK_C8FCBCFB_8280_4F1D_9159_8E1C3F95E045

mLatSubset *qlua_checkLatSubset(lua_State *L, int idx, mLattice *S);
mLatSubset *qlua_newLatSubset(lua_State *L, int Sidx);

int init_latsubset(lua_State *L);
int fini_latsubset(lua_State *L);

#endif /* !defined(MARK_C8FCBCFB_8280_4F1D_9159_8E1C3F95E045) */
