#ifndef MARK_DB3FC6CE_E1BB_4A19_924C_FD12859F06F9
#define MARK_DB3FC6CE_E1BB_4A19_924C_FD12859F06F9

#if USE_Nc2
typedef struct {
    void *ptr;
} mSeqColVec2;

mSeqColVec2 *qlua_checkSeqColVec2(lua_State *L, int idx, int nc);
mSeqColVec2 *qlua_newSeqColVec2(lua_State *L, int nc);
mSeqColVec2 *qlua_newZeroSeqColVec2(lua_State *L, int nc);
#endif

#if USE_Nc3
typedef struct {
    void *ptr;
} mSeqColVec3;

mSeqColVec3 *qlua_checkSeqColVec3(lua_State *L, int idx, int nc);
mSeqColVec3 *qlua_newSeqColVec3(lua_State *L, int nc);
mSeqColVec3 *qlua_newZeroSeqColVec3(lua_State *L, int nc);
#endif

#if USE_NcN
typedef struct {
    int nc;
    void *ptr;
} mSeqColVecN;

mSeqColVecN *qlua_checkSeqColVecN(lua_State *L, int idx, int nc);
mSeqColVecN *qlua_newSeqColVecN(lua_State *L, int nc);
mSeqColVecN *qlua_newZeroSeqColVecN(lua_State *L, int nc);
#endif

int init_seqcolvec(lua_State *L);
void fini_seqcolvec(void);

int q_v_gaussian_N(lua_State *L);

#endif /* !defined(MARK_DB3FC6CE_E1BB_4A19_924C_FD12859F06F9) */
