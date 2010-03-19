#ifndef MARK_02FBE25A_A3B1_412A_9216_C523C560ED1C
#define MARK_02FBE25A_A3B1_412A_9216_C523C560ED1C

typedef struct {
    int       size;
    QLA_Int  *idx;
} mLatMulti;

extern const char mtnLatMultSet[];

mLatMulti *qlua_checkLatMulti(lua_State *L, int idx, mLattice *S);

int init_latmulti(lua_State *L);
int fini_latmulti(lua_State *L);

#endif /* !defined(MARK_02FBE25A_A3B1_412A_9216_C523C560ED1C) */
