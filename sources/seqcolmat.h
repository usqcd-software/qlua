#ifndef MARK_3EA0F7A6_DE4E_422E_BF5A_73999D4AFEB7
#define MARK_3EA0F7A6_DE4E_422E_BF5A_73999D4AFEB7

#if USE_Nc2
typedef struct {
    void *ptr;
} mSeqColMat2;

mSeqColMat2 *qlua_checkSeqColMat2(lua_State *L, int idx, int nc);
mSeqColMat2 *qlua_newSeqColMat2(lua_State *L, int nc);
mSeqColMat2 *qlua_newZeroSeqColMat2(lua_State *L, int nc);
#endif

#if USE_Nc3
typedef struct {
    void *ptr;
} mSeqColMat3;

mSeqColMat3 *qlua_checkSeqColMat3(lua_State *L, int idx, int nc);
mSeqColMat3 *qlua_newSeqColMat3(lua_State *L, int nc);
mSeqColMat3 *qlua_newZeroSeqColMat3(lua_State *L, int nc);
#endif

#if USE_NcN
typedef struct {
    int nc;
    void *ptr;
} mSeqColMatN;

mSeqColMatN *qlua_checkSeqColMatN(lua_State *L, int idx, int nc);
mSeqColMatN *qlua_newSeqColMatN(lua_State *L, int nc);
mSeqColMatN *qlua_newZeroSeqColMatN(lua_State *L, int nc);
#endif

int init_seqcolmat(lua_State *L);
int fini_seqcolmat(lua_State *L);

int q_m_gaussian_N(lua_State *L);

#endif /* !defined(MARK_3EA0F7A6_DE4E_422E_BF5A_73999D4AFEB7) */
