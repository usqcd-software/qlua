#ifndef MARK_B05C8FF9_2888_4CF2_B9AD_29743658DC2A
#define MARK_B05C8FF9_2888_4CF2_B9AD_29743658DC2A

typedef struct {
    QLA_RandomState state;
} mSeqRandom;

mSeqRandom *qlua_checkSeqRandom(lua_State *L, int idx);
mSeqRandom *qlua_newSeqRandom(lua_State *L);

int init_seqrandom(lua_State *L);
void fini_seqrandom(void);

#endif /* !defined(MARK_B05C8FF9_2888_4CF2_B9AD_29743658DC2A) */
