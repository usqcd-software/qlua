#ifndef MARK_73980191_6DDF_41BA_ADD9_143100E683FF
#define MARK_73980191_6DDF_41BA_ADD9_143100E683FF

typedef struct {
    int size;
    QLA_Int val[1];
} mVecInt;

typedef struct {
    int size;
    QLA_D_Real val[1];
} mVecReal;

typedef struct {
    int size;
    QLA_D_Complex val[1];
} mVecComplex;

extern const char mtnVecInt[];
extern const char mtnVecReal[];
extern const char mtnVecComplex[];

int init_vector(lua_State *L);
int fini_vector(lua_State *L);

mVecInt *qlua_checkVecInt(lua_State *L, int idx);
mVecInt *qlua_newVecInt(lua_State *L, int size);
mVecReal *qlua_checkVecReal(lua_State *L, int idx) ;
mVecReal *qlua_newVecReal(lua_State *L, int size);
mVecComplex *qlua_checkVecComplex(lua_State *L, int idx);
mVecComplex *qlua_newVecComplex(lua_State *L, int size);

#endif /* !defined(MARK_73980191_6DDF_41BA_ADD9_143100E683FF) */
