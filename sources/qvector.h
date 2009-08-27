#ifndef MARK_73980191_6DDF_41BA_ADD9_143100E683FF
#define MARK_73980191_6DDF_41BA_ADD9_143100E683FF

typedef struct {
    int size;
    QLA_Int val[1];
} tVecInt;

typedef struct {
    int size;
    QLA_D_Real val[1];
} tVecDouble;

typedef struct {
    int size;
    QLA_D_Complex val[1];
} tVecComplex;

extern const char *mtnVecInt;
extern const char *mtnVecDouble;
extern const char *mtnVecComplex;

int init_vector(lua_State *L);
int fini_vector(lua_State *L);

tVecInt *q_checkVecInt(lua_State *L, int idx);
tVecInt *q_newVecInt(lua_State *L, int size);
tVecDouble *q_checkVecDouble(lua_State *L, int idx) ;
tVecDouble *q_newVecDouble(lua_State *L, int size);
tVecComplex *q_checkVecComplex(lua_State *L, int idx);
tVecComplex *q_newVecComplex(lua_State *L, int size);

#endif /* !defined(MARK_73980191_6DDF_41BA_ADD9_143100E683FF) */
