#ifndef MARK_745BEFA9_4A4B_4C2F_AF23_62D5F3F05CDF
#define MARK_745BEFA9_4A4B_4C2F_AF23_62D5F3F05CDF

#if USE_Nc2
typedef struct {
    void *ptr;
} mSeqDirFerm2;

mSeqDirFerm2 *qlua_checkSeqDirFerm2(lua_State *L, int idx, int nc);
mSeqDirFerm2 *qlua_newSeqDirFerm2(lua_State *L, int nc);
mSeqDirFerm2 *qlua_newZeroSeqDirFerm2(lua_State *L, int nc);
#endif

#if USE_Nc3
typedef struct {
    void *ptr;
} mSeqDirFerm3;

mSeqDirFerm3 *qlua_checkSeqDirFerm3(lua_State *L, int idx, int nc);
mSeqDirFerm3 *qlua_newSeqDirFerm3(lua_State *L, int nc);
mSeqDirFerm3 *qlua_newZeroSeqDirFerm3(lua_State *L, int nc);
#endif

#if USE_NcN
typedef struct {
    int nc;
    void *ptr;
} mSeqDirFermN;

mSeqDirFermN *qlua_checkSeqDirFermN(lua_State *L, int idx, int nc);
mSeqDirFermN *qlua_newSeqDirFermN(lua_State *L, int nc);
mSeqDirFermN *qlua_newZeroSeqDirFermN(lua_State *L, int nc);
#endif

int init_seqdirferm(lua_State *L);
void fini_seqdirferm(void);

int q_d_gaussian_N(lua_State *L);

#endif /* !defined(MARK_745BEFA9_4A4B_4C2F_AF23_62D5F3F05CDF) */
