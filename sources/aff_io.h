#ifndef MARK_B3F5C1AE_EF2E_4D60_89EF_56F7D67F4287
#define MARK_B3F5C1AE_EF2E_4D60_89EF_56F7D67F4287

typedef struct {
    struct AffReader_s  *ptr;
    struct AffNode_s    *dir;
} mAffReader;

typedef struct {
    struct AffWriter_s  *ptr;
    struct AffNode_s    *dir;
} mAffWriter;

int init_aff_io(lua_State *L);
int fini_aff_io(lua_State *L);

mAffReader *qlua_checkAffReader(lua_State *L, int idx);
mAffWriter *qlua_checkAffWriter(lua_State *L, int idx);

#endif /* !defined(MARK_B3F5C1AE_EF2E_4D60_89EF_56F7D67F4287) */
