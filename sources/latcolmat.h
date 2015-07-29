#ifndef MARK_E53A3544_2065_4CA4_888F_72B524D3269B
#define MARK_E53A3544_2065_4CA4_888F_72B524D3269B

#if USE_Nc2
typedef struct {
    QDP_D2_ColorMatrix *ptr;
} mLatColMat2;

mLatColMat2 *qlua_checkLatColMat2(lua_State *L, int idx, mLattice *S, int nc);
mLatColMat2 *qlua_newLatColMat2(lua_State *L, int Sidx, int nc);
mLatColMat2 *qlua_newZeroLatColMat2(lua_State *L, int Sidx, int nc);
#endif

#if USE_Nc3
typedef struct {
    QDP_D3_ColorMatrix *ptr;
} mLatColMat3;

mLatColMat3 *qlua_checkLatColMat3(lua_State *L, int idx, mLattice *S, int nc);
mLatColMat3 *qlua_newLatColMat3(lua_State *L, int Sidx, int nc);
mLatColMat3 *qlua_newZeroLatColMat3(lua_State *L, int Sidx, int nc);
#endif

#if USE_NcN
typedef struct {
    int nc;
    QDP_DN_ColorMatrix *ptr;
} mLatColMatN;

mLatColMatN *qlua_checkLatColMatN(lua_State *L, int idx, mLattice *S, int nc);
mLatColMatN *qlua_newLatColMatN(lua_State *L, int Sidx, int nc);
mLatColMatN *qlua_newZeroLatColMatN(lua_State *L, int Sidx, int nc);
#endif

int init_latcolmat(lua_State *L);
void fini_latcolmat(void);

int q_M_gaussian(lua_State *L);
int q_M_gaussian_N(lua_State *L);

#endif /* !defined(MARK_E53A3544_2065_4CA4_888F_72B524D3269B) */
