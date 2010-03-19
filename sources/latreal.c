#include "qlua.h"                                                    /* DEPS */
#include "qvector.h"                                                 /* DEPS */
#include "latsubset.h"                                               /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "latmulti.h"                                                /* DEPS */
#include "latreal.h"                                                 /* DEPS */
#include "latint.h"                                                  /* DEPS */
#include "latrandom.h"                                               /* DEPS */
#include "latcomplex.h"                                              /* DEPS */
#include "qmp.h"
#include "qdp_df.h"

static int
q_RD_add_RF(lua_State *L)
{
    mLatRealD *a = qlua_checkLatRealD(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatRealF *b = qlua_checkLatRealF(L, 2, S);
    mLatRealD *c = qlua_newLatRealD(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_DF_R_eq_R(c->ptr, b->ptr, *S->qss);
    QDP_D_R_peq_R(c->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_RF_add_RD(lua_State *L)
{
    mLatRealF *a = qlua_checkLatRealF(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatRealD *b = qlua_checkLatRealD(L, 2, S);
    mLatRealD *c = qlua_newLatRealD(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_DF_R_eq_R(c->ptr, a->ptr, *S->qss);
    QDP_D_R_peq_R(c->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_RD_sub_RF(lua_State *L)
{
    mLatRealD *a = qlua_checkLatRealD(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatRealF *b = qlua_checkLatRealF(L, 2, S);
    int Sidx = lua_gettop(L);
    mLatRealD *x = qlua_newLatRealD(L, Sidx);
    mLatRealD *c = qlua_newLatRealD(L, Sidx);

    CALL_QDP(L);
    QDP_DF_R_eq_R(x->ptr, b->ptr, *S->qss);
    QDP_D_R_eq_R_minus_R(c->ptr, a->ptr, x->ptr, *S->qss);

    return 1;
}

static int
q_RF_sub_RD(lua_State *L)
{
    mLatRealF *a = qlua_checkLatRealF(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatRealD *b = qlua_checkLatRealD(L, 2, S);
    mLatRealD *c = qlua_newLatRealD(L, lua_gettop(L));

    CALL_QDP(L);
    QDP_DF_R_eq_R(c->ptr, a->ptr, *S->qss);
    QDP_D_R_meq_R(c->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_RD_mul_RF(lua_State *L)
{
    mLatRealD *a = qlua_checkLatRealD(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatRealF *b = qlua_checkLatRealF(L, 2, S);
    int Sidx = lua_gettop(L);
    mLatRealD *x = qlua_newLatRealD(L, Sidx);
    mLatRealD *c = qlua_newLatRealD(L, Sidx);

    CALL_QDP(L);
    QDP_DF_R_eq_R(x->ptr, b->ptr, *S->qss);
    QDP_D_R_eq_R_times_R(c->ptr, a->ptr, x->ptr, *S->qss);

    return 1;
}

static int
q_RF_mul_RD(lua_State *L)
{
    mLatRealF *a = qlua_checkLatRealF(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatRealD *b = qlua_checkLatRealD(L, 2, S);
    int Sidx = lua_gettop(L);
    mLatRealD *x = qlua_newLatRealD(L, Sidx);
    mLatRealD *c = qlua_newLatRealD(L, Sidx);

    CALL_QDP(L);
    QDP_DF_R_eq_R(x->ptr, a->ptr, *S->qss);
    QDP_D_R_eq_R_times_R(c->ptr, x->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_RD_div_RF(lua_State *L)
{
    mLatRealD *a = qlua_checkLatRealD(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatRealF *b = qlua_checkLatRealF(L, 2, S);
    int Sidx = lua_gettop(L);
    mLatRealD *x = qlua_newLatRealD(L, Sidx);
    mLatRealD *c = qlua_newLatRealD(L, Sidx);

    CALL_QDP(L);
    QDP_DF_R_eq_R(x->ptr, b->ptr, *S->qss);
    QDP_D_R_eq_R_divide_R(c->ptr, a->ptr, x->ptr, *S->qss);

    return 1;
}

static int
q_RF_div_RD(lua_State *L)
{
    mLatRealF *a = qlua_checkLatRealF(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatRealD *b = qlua_checkLatRealD(L, 2, S);
    int Sidx = lua_gettop(L);
    mLatRealD *x = qlua_newLatRealD(L, Sidx);
    mLatRealD *c = qlua_newLatRealD(L, Sidx);

    CALL_QDP(L);
    QDP_DF_R_eq_R(x->ptr, a->ptr, *S->qss);
    QDP_D_R_eq_R_divide_R(c->ptr, x->ptr, b->ptr, *S->qss);

    return 1;
}

static int
q_RD_float(lua_State *L)
{
    mLatRealD *a = qlua_checkLatRealD(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatRealF *c = qlua_newLatRealF(L, lua_gettop(L));
    
    CALL_QDP(L);
    QDP_FD_R_eq_R(c->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_RD_double(lua_State *L)
{
    mLatRealD *a = qlua_checkLatRealD(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatRealD *c = qlua_newLatRealD(L, lua_gettop(L));
    
    CALL_QDP(L);
    QDP_D_R_eq_R(c->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_RF_float(lua_State *L)
{
    mLatRealF *a = qlua_checkLatRealF(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatRealF *c = qlua_newLatRealF(L, lua_gettop(L));
    
    CALL_QDP(L);
    QDP_F_R_eq_R(c->ptr, a->ptr, *S->qss);

    return 1;
}

static int
q_RF_double(lua_State *L)
{
    mLatRealF *a = qlua_checkLatRealF(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    mLatRealD *c = qlua_newLatRealD(L, lua_gettop(L));
    
    CALL_QDP(L);
    QDP_DF_R_eq_R(c->ptr, a->ptr, *S->qss);

    return 1;
}

#define Qab(p,s)   p ## D ## s ## D
#define Qa1(p)     p ## D
#define Qa2(p,s)   p ## D ## s
#define Qx         'D'
#include "latreal-x.c"                                              /* DEPS */

#define Qab(p,s)   p ## F ## s ## F
#define Qa1(p)     p ## F
#define Qa2(p,s)   p ## F ## s
#define Qx         'F'
#include "latreal-x.c"                                              /* DEPS */

static int
q_latreal(lua_State *L)
{
    mLattice *S = qlua_checkLattice(L, 1);

    switch (S->precision) {
    case 'F':
        return q_latrealF(L);
    case 'D':
        return q_latrealD(L);
    default:
        return luaL_error(L, "*** internal error: unsupported precision ***");
    }
}

int
q_R_random(lua_State *L)
{
    mLattice *S = qlua_ObjLattice(L, 1);

    switch (S->precision) {
    case 'F':
        return q_RF_random(L);
    case 'D':
        return q_RD_random(L);
    default:
        return luaL_error(L, "*** internal error: unsupported precision ***");
    }
}

int
q_R_gaussian(lua_State *L)
{
    mLattice *S = qlua_ObjLattice(L, 1);

    switch (S->precision) {
    case 'F':
        return q_RF_gaussian(L);
    case 'D':
        return q_RD_gaussian(L);
    default:
        return luaL_error(L, "*** internal error: unsupported precision ***");
    }
}


static struct luaL_Reg fLatReal[] = {
    { "RealD",      q_latrealD },
    { "RealF",      q_latrealF },
    { "Real",       q_latreal },
    { NULL,         NULL }
};

int
init_latreal(lua_State *L)
{
    static const QLUA_Op2 ops[] = {
        {qlua_add_table, qLatRealD, qLatRealD, q_RD_add_RD },
        {qlua_add_table, qLatRealD, qLatRealF, q_RD_add_RF },
        {qlua_add_table, qLatRealF, qLatRealF, q_RF_add_RF },
        {qlua_add_table, qLatRealF, qLatRealD, q_RF_add_RD },
        {qlua_add_table, qReal,     qLatRealD, q_r_add_RD  },
        {qlua_add_table, qReal,     qLatRealF, q_r_add_RF  },
        {qlua_add_table, qLatRealD, qReal,     q_RD_add_r  },
        {qlua_add_table, qLatRealF, qReal,     q_RF_add_r  },
        {qlua_sub_table, qLatRealD, qLatRealD, q_RD_sub_RD },
        {qlua_sub_table, qLatRealD, qLatRealF, q_RD_sub_RF },
        {qlua_sub_table, qLatRealF, qLatRealF, q_RF_sub_RF },
        {qlua_sub_table, qLatRealF, qLatRealD, q_RF_sub_RD },
        {qlua_sub_table, qReal,     qLatRealD, q_r_sub_RD  },
        {qlua_sub_table, qReal,     qLatRealF, q_r_sub_RF  },
        {qlua_sub_table, qLatRealD, qReal,     q_RD_sub_r  },
        {qlua_sub_table, qLatRealF, qReal,     q_RF_sub_r  },
        {qlua_mul_table, qLatRealD, qLatRealD, q_RD_mul_RD },
        {qlua_mul_table, qLatRealD, qLatRealF, q_RD_mul_RF },
        {qlua_mul_table, qLatRealF, qLatRealF, q_RF_mul_RF },
        {qlua_mul_table, qLatRealF, qLatRealD, q_RF_mul_RD },
        {qlua_mul_table, qReal,     qLatRealD, q_r_mul_RD  },
        {qlua_mul_table, qReal,     qLatRealF, q_r_mul_RF  },
        {qlua_mul_table, qLatRealD, qReal,     q_RD_mul_r  },
        {qlua_mul_table, qLatRealF, qReal,     q_RF_mul_r  },
        {qlua_div_table, qLatRealD, qLatRealD, q_RD_div_RD },
        {qlua_div_table, qLatRealD, qLatRealF, q_RD_div_RF },
        {qlua_div_table, qLatRealF, qLatRealF, q_RF_div_RF },
        {qlua_div_table, qLatRealF, qLatRealD, q_RF_div_RD },
        {qlua_div_table, qReal,     qLatRealD, q_r_div_RD  },
        {qlua_div_table, qReal,     qLatRealF, q_r_div_RF  },
        {qlua_div_table, qLatRealD, qReal,     q_RD_div_r  },
        {qlua_div_table, qLatRealF, qReal,     q_RF_div_r  },
        {NULL,           qNoType,   qNoType,   NULL        }
    };

    luaL_getmetatable(L, opLattice);
    luaL_register(L, NULL, fLatReal);
    lua_pop(L, 1);
    qlua_reg_op2(ops);
    qlua_reg_dot(qLatRealD, q_RD_mul_RD);
    qlua_reg_dot(qLatRealF, q_RF_mul_RF);

    return 0;
}

int
fini_latreal(lua_State *L)
{
    return 0;
}
