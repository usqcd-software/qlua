#ifndef MARK_1EF46601_F4FC_44FF_8701_E58341017024
#define MARK_1EF46601_F4FC_44FF_8701_E58341017024

typedef struct {
    QDP_F_Complex *ptr;
} mLatComplexF;

typedef struct {
    QDP_D_Complex *ptr;
} mLatComplexD;

int init_latcomplex(lua_State *L);
int fini_latcomplex(lua_State *L);

mLatComplexD *qlua_checkLatComplexD(lua_State *L, int idx, mLattice *S);
mLatComplexF *qlua_checkLatComplexF(lua_State *L, int idx, mLattice *S);
mLatComplexD *qlua_newLatComplexD(lua_State *L, int S_idx);
mLatComplexF *qlua_newLatComplexF(lua_State *L, int S_idx);

int q_C_gaussian(lua_State *L);

#endif /* !defined(MARK_1EF46601_F4FC_44FF_8701_E58341017024) */
