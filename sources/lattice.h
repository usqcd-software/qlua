#ifndef MARK_B3A3A85C_2C48_4B4C_B926_440CBD0CB411
#define MARK_B3A3A85C_2C48_4B4C_B926_440CBD0CB411

extern int qRank;
extern int *qDim;

int init_lattice(lua_State *L);
int fini_lattice(lua_State *L);

int *qlua_latcoord(lua_State *L, int idx);                   /* lc[] or NULL */
int *qlua_checklatcoord(lua_State *L, int idx);

#endif /* !defined(MARK_B3A3A85C_2C48_4B4C_B926_440CBD0CB411) */
