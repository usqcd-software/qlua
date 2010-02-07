#ifndef MARK_FE647744_2081_4E52_9AE1_899F1FA89AE5
#define MARK_FE647744_2081_4E52_9AE1_899F1FA89AE5

typedef struct {
    int size_l;
    int size_r;
    double val[1];
} mMatReal;

typedef struct {
    int size_l;
    int size_r;
    double val[2];
} mMatComplex;

extern const char mtnMatReal[];
extern const char mtnMatComplex[];

int init_matrix(lua_State *L);
int fini_matrix(lua_State *L);

mMatReal *qlua_checkMatReal(lua_State *L, int idx) ;
mMatReal *qlua_newMatReal(lua_State *L, int size_l, int size_r);
mMatComplex *qlua_checkMatComplex(lua_State *L, int idx) ;
mMatComplex *qlua_newMatComplex(lua_State *L, int size_l, int size_r);

#endif /* !defined(MARK_FE647744_2081_4E52_9AE1_899F1FA89AE5) */
