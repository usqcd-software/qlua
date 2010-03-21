#ifndef MARK_2D398035_F07C_4356_BAC3_5795771D71FD
#define MARK_2D398035_F07C_4356_BAC3_5795771D71FD

typedef struct {
    QDP_D2_ColorVector *ptr;
} mLatColVec2;

typedef struct {
    QDP_D3_ColorVector *ptr;
} mLatColVec3;

typedef struct {
    int nc;
    QDP_DN_ColorVector *ptr;
} mLatColVecN;

int init_latcolvec(lua_State *L);
int fini_latcolvec(lua_State *L);

mLatColVec2 *qlua_checkLatColVec2(lua_State *L, int idx, mLattice *S, int nc);
mLatColVec3 *qlua_checkLatColVec3(lua_State *L, int idx, mLattice *S, int nc);
mLatColVecN *qlua_checkLatColVecN(lua_State *L, int idx, mLattice *S, int nc);
mLatColVec2 *qlua_newLatColVec2(lua_State *L, int Sidx, int nc);
mLatColVec3 *qlua_newLatColVec3(lua_State *L, int Sidx, int nc);
mLatColVecN *qlua_newLatColVecN(lua_State *L, int Sidx, int nc);

int q_V_gaussian(lua_State *L);
int q_V_gaussian_N(lua_State *L);

#endif /* !defined(MARK_2D398035_F07C_4356_BAC3_5795771D71FD) */
