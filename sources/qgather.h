#ifndef MARK_82886201_7371_4C6B_B39F_29E0C238A18A
#define MARK_82886201_7371_4C6B_B39F_29E0C238A18A

typedef struct {
    QDP_Collect *ptr; /* we had to use collect/scatter instead of gather/scatter in QDP ... */
} mGather;

int init_gather(lua_State *L);
void fini_gather(void);

mGather *qlua_checkGather(lua_State *L, int idx);

#endif /* !defined(MARK_82886201_7371_4C6B_B39F_29E0C238A18A) */
