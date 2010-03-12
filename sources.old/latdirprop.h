#ifndef MARK_5F7EAA43_B7D7_4A36_B574_D9707797522B
#define MARK_5F7EAA43_B7D7_4A36_B574_D9707797522B

typedef struct {
    QDP_DiracPropagator *ptr;
} mLatDirProp;

extern const char mtnLatDirProp[];

int init_latdirprop(lua_State *L);
int fini_latdirprop(lua_State *L);

mLatDirProp *qlua_checkLatDirProp(lua_State *L, int idx);
mLatDirProp *qlua_newLatDirProp(lua_State *L);

int q_P_gaussian(lua_State *L);

#endif /* !defined(MARK_5F7EAA43_B7D7_4A36_B574_D9707797522B) */
