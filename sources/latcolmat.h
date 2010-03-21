#ifndef MARK_E53A3544_2065_4CA4_888F_72B524D3269B
#define MARK_E53A3544_2065_4CA4_888F_72B524D3269B

typedef struct {
    QDP_D2_ColorMatrix *ptr;
} mLatColMat2;

typedef struct {
    QDP_D3_ColorMatrix *ptr;
} mLatColMat3;

typedef struct {
    int nc;
    QDP_DN_ColorMatrix *ptr;
} mLatColMatN;

int init_latcolmat(lua_State *L);
int fini_latcolmat(lua_State *L);

mLatColMat2 *qlua_checkLatColMat2(lua_State *L, int idx, mLattice *S, int nc);
mLatColMat3 *qlua_checkLatColMat3(lua_State *L, int idx, mLattice *S, int nc);
mLatColMatN *qlua_checkLatColMatN(lua_State *L, int idx, mLattice *S, int nc);
mLatColMat2 *qlua_newLatColMat2(lua_State *L, int Sidx, int nc);
mLatColMat3 *qlua_newLatColMat3(lua_State *L, int Sidx, int nc);
mLatColMatN *qlua_newLatColMatN(lua_State *L, int Sidx, int nc);

int q_M_gaussian(lua_State *L);
int q_M_gaussian_N(lua_State *L);

#endif /* !defined(MARK_E53A3544_2065_4CA4_888F_72B524D3269B) */
