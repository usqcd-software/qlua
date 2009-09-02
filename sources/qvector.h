#ifndef MARK_73980191_6DDF_41BA_ADD9_143100E683FF
#define MARK_73980191_6DDF_41BA_ADD9_143100E683FF

typedef struct {
    int size;
    QLA_Int val[1];
} tVecInt;

typedef struct {
    int size;
    QLA_D_Real val[1];
} tVecReal;

typedef struct {
    int size;
    QLA_D_Complex val[1];
} tVecComplex;

extern const char mtnVecInt[];
extern const char mtnVecReal[];
extern const char mtnVecComplex[];

int init_vector(lua_State *L);
int fini_vector(lua_State *L);

tVecInt *qlua_checkVecInt(lua_State *L, int idx);
tVecInt *qlua_newVecInt(lua_State *L, int size);
tVecReal *qlua_checkVecReal(lua_State *L, int idx) ;
tVecReal *qlua_newVecReal(lua_State *L, int size);
tVecComplex *qlua_checkVecComplex(lua_State *L, int idx);
tVecComplex *qlua_newVecComplex(lua_State *L, int size);

#endif /* !defined(MARK_73980191_6DDF_41BA_ADD9_143100E683FF) */
