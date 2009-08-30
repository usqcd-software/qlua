#ifndef MARK_2D398035_F07C_4356_BAC3_5795771D71FD
#define MARK_2D398035_F07C_4356_BAC3_5795771D71FD

typedef struct {
    QDP_ColorVector *ptr;
} mLatColVec;

extern const char mtnLatColVec[];

int init_latcolvec(lua_State *L);
int fini_latcolvec(lua_State *L);

mLatColVec *qlua_checkLatColVec(lua_State *L, int idx);
mLatColVec *qlua_newLatColVec(lua_State *L);

int q_V_gaussian(lua_State *L);

int q_V_dot(lua_State *L);
int q_V_add_V(lua_State *L);
int q_V_sub_V(lua_State *L);
int q_r_mul_V(lua_State *L);
int q_V_mul_r(lua_State *L);
int q_c_mul_V(lua_State *L);
int q_V_mul_c(lua_State *L);


#endif /* !defined(MARK_2D398035_F07C_4356_BAC3_5795771D71FD) */
