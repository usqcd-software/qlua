#ifndef MARK_B3A3A85C_2C48_4B4C_B926_440CBD0CB411
#define MARK_B3A3A85C_2C48_4B4C_B926_440CBD0CB411

extern const char opLattice[];
extern int qRank;
extern int *qDim;
extern QDP_Subset qCurrent;

int *qlua_latcoord(lua_State *L, int idx);                   /* lc[] or NULL */
int *qlua_checklatcoord(lua_State *L, int idx);

QDP_Shift qlua_checkShift(lua_State *L, int idx);
QDP_ShiftDir qlua_checkShiftDir(lua_State *L, int idx);

int init_lattice(lua_State *L);
int fini_lattice(lua_State *L);

#endif /* !defined(MARK_B3A3A85C_2C48_4B4C_B926_440CBD0CB411) */
