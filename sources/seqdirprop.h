#ifndef MARK_EB604AB4_371D_46F4_94DD_83057F574BAC
#define MARK_EB604AB4_371D_46F4_94DD_83057F574BAC

#if USE_Nc2
typedef struct {
    void *ptr;
} mSeqDirProp2;

mSeqDirProp2 *qlua_checkSeqDirProp2(lua_State *L, int idx, int nc);
mSeqDirProp2 *qlua_newSeqDirProp2(lua_State *L, int nc);
mSeqDirProp2 *qlua_newZeroSeqDirProp2(lua_State *L, int nc);
#endif

#if USE_Nc3
typedef struct {
    void *ptr;
} mSeqDirProp3;

mSeqDirProp3 *qlua_checkSeqDirProp3(lua_State *L, int idx, int nc);
mSeqDirProp3 *qlua_newSeqDirProp3(lua_State *L, int nc);
mSeqDirProp3 *qlua_newZeroSeqDirProp3(lua_State *L, int nc);
#endif

#if USE_NcN
typedef struct {
    int nc;
    void *ptr;
} mSeqDirPropN;

mSeqDirPropN *qlua_checkSeqDirPropN(lua_State *L, int idx, int nc);
mSeqDirPropN *qlua_newSeqDirPropN(lua_State *L, int nc);
mSeqDirPropN *qlua_newZeroSeqDirPropN(lua_State *L, int nc);
#endif

int init_seqdirprop(lua_State *L);
void fini_seqdirprop(void);

int q_p_gaussian_N(lua_State *L);

#endif /* !defined(MARK_EB604AB4_371D_46F4_94DD_83057F574BAC) */
