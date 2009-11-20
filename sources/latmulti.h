#ifndef MARK_02FBE25A_A3B1_412A_9216_C523C560ED1C
#define MARK_02FBE25A_A3B1_412A_9216_C523C560ED1C

typedef struct {
    int         axis;
    int         count;
    void       *arg;
    QDP_Subset *subset;
} mLatMulti;

extern const char mtnLatMultSet[];

void qlua_checkLatMulti(lua_State *L, int idx);
int qlua_LatMultiSize(lua_State *L, int idx);
mLatInt *qlua_LatMultiIndex(lua_State *L, int idx); /* see source */

int init_latmulti(lua_State *L);
int fini_latmulti(lua_State *L);

#endif /* !defined(MARK_02FBE25A_A3B1_412A_9216_C523C560ED1C) */
