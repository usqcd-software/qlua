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
    QDP_Lattice *lat;
    int net_forced;
    lua_State *L;
    int id;
    int nc;
    int rank;
    int *dim;
    int node;
    int *net;
    int *neighbor_up;
    int *neighbor_down;
    QDP_Subset *qss;
    mLatSubset lss;
    QDP_Subset *none;
    QDP_Subset all;
    QDP_Subset even;
    QDP_Subset odd;
} mLattice;

extern const char opLattice[];

int *qlua_latcoord(lua_State *L, int idx, mLattice *S);      /* lc[] or NULL */
int *qlua_checklatcoord(lua_State *L, int idx, mLattice *S);
void qlua_verifylatcoord(lua_State *L, int *coord, mLattice *S);

mLattice *qlua_checkLattice(lua_State *L, int idx);
mLattice *qlua_ObjLattice(lua_State *L, int idx);

QDP_Shift qlua_checkShift(lua_State *L, int idx, mLattice *S);
QDP_ShiftDir qlua_checkShiftDir(lua_State *L, int idx);

int init_lattice(lua_State *L);
void fini_lattice(void);

#endif /* !defined(MARK_B3A3A85C_2C48_4B4C_B926_440CBD0CB411) */
