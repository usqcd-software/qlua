#ifndef MARK_E53A3544_2065_4CA4_888F_72B524D3269B
#define MARK_E53A3544_2065_4CA4_888F_72B524D3269B

typedef struct {
    QDP_ColorMatrix *ptr;
} mLatColMat;

extern const char mtnLatColMat[];

int init_latcolmat(lua_State *L);
int fini_latcolmat(lua_State *L);

mLatColMat *qlua_checkLatColMat(lua_State *L, int idx);
mLatColMat *qlua_newLatColMat(lua_State *L);

int q_M_gaussian(lua_State *L);

#endif /* !defined(MARK_E53A3544_2065_4CA4_888F_72B524D3269B) */
