#ifndef MARK_40238D8A_7D6D_47BA_9C1B_61A31864BFFA
#define MARK_40238D8A_7D6D_47BA_9C1B_61A31864BFFA

typedef struct {
    QDP_Scatter *ptr;
} mScatter;

int init_scatter(lua_State *L);
void fini_scatter(void);

mScatter *qlua_checkScatter(lua_State *L, int idx);

#endif /* !defined(MARK_40238D8A_7D6D_47BA_9C1B_61A31864BFFA) */
