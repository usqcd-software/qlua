#ifndef MARK_1EF46601_F4FC_44FF_8701_E58341017024
#define MARK_1EF46601_F4FC_44FF_8701_E58341017024

typedef struct {
    QDP_Complex *ptr;
} mLatComplex;

int init_latcomplex(lua_State *L);
void fini_latcomplex(void);

mLatComplex *qlua_checkLatComplex(lua_State *L, int idx, mLattice *S);
mLatComplex *qlua_newLatComplex(lua_State *L, int S_idx);
mLatComplex *qlua_newZeroLatComplex(lua_State *L, int S_idx);

int q_C_gaussian(lua_State *L);

#endif /* !defined(MARK_1EF46601_F4FC_44FF_8701_E58341017024) */
