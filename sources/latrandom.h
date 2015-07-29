#ifndef MARK_3DF9923F_966E_47B3_ABAA_79CFF3A58E8C
#define MARK_3DF9923F_966E_47B3_ABAA_79CFF3A58E8C

typedef struct {
    QDP_RandomState *ptr;
} mLatRandom;

mLatRandom *qlua_checkLatRandom(lua_State *L, int idx, mLattice *S);
mLatRandom *qlua_newLatRandom(lua_State *L, int Sidx);
mLatRandom *qlua_newZeroLatRandom(lua_State *L, int Sidx);

int init_latrandom(lua_State *L);
void fini_latrandom(void);

#endif /* !defined(MARK_3DF9923F_966E_47B3_ABAA_79CFF3A58E8C) */
