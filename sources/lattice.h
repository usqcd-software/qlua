#ifndef MARK_B3A3A85C_2C48_4B4C_B926_440CBD0CB411
#define MARK_B3A3A85C_2C48_4B4C_B926_440CBD0CB411

typedef enum {
    qss_all,
    qss_even,
    qss_odd,
    qss_none,
    qss_last_static = qss_none,
    qss_slice,
    qss_upper,
    qss_lower
} qSubsetClass;

typedef struct {
    qSubsetClass  cl;
    int           axis;
    int           position;
    QDP_Int      *mask;
} mLatSubset;

typedef struct {
    int id;
    int nc;
    int rank;
    int *dim;
    QDP_Subset *qss;
    mLatSubset lss;
} mLattice;

extern const char opLattice[];

int *qlua_latcoord(lua_State *L, int idx, mLattice *S);      /* lc[] or NULL */
int *qlua_checklatcoord(lua_State *L, int idx, mLattice *S);
void qlua_verifylatcoord(lua_State *L, int *coord, mLattice *S);
int *qlua_intarray(lua_State *L, int idx, int *out_dim);
int *qlua_checkintarray(lua_State *L, int idx, int dim, int *out_dim);

mLattice *qlua_checkLattice(lua_State *L, int idx);
mLattice *qlua_ObjLattice(lua_State *L, int idx);

QDP_Shift qlua_checkShift(lua_State *L, int idx, mLattice *S);
QDP_ShiftDir qlua_checkShiftDir(lua_State *L, int idx);

int init_lattice(lua_State *L);
int fini_lattice(lua_State *L);

#endif /* !defined(MARK_B3A3A85C_2C48_4B4C_B926_440CBD0CB411) */
