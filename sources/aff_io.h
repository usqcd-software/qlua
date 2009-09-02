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

void qlua_Aff_enter(lua_State *L);
void qlua_Aff_leave(void);
struct AffNode_s *qlua_AffReaderChPath(mAffReader *b, const char *p);
struct AffNode_s *qlua_AffWriterMkPath(mAffWriter *b, const char *p);

#endif /* !defined(MARK_B3F5C1AE_EF2E_4D60_89EF_56F7D67F4287) */
