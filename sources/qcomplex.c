#include "qlua.h"                                                    /* DEPS */
#include "qcomplex.h"                                                /* DEPS */
#include "seqrandom.h"                                               /* DEPS */

const char mtnComplex[] = "qlua.mtComplex";

QLA_D_Complex *
qlua_checkComplex(lua_State *L, int idx)
{
    void *ud = luaL_checkudata(L, idx, mtnComplex);

    luaL_argcheck(L, ud != NULL, idx, "complex expected");

    return (QLA_D_Complex *)ud;
}

QLA_D_Complex *
qlua_newComplex(lua_State *L)
{
    QLA_D_Complex *z = (QLA_D_Complex *)lua_newuserdata(L, sizeof (QLA_D_Complex));
    
    luaL_getmetatable(L, mtnComplex);
    lua_setmetatable(L, -2);

    return z;
}


static int
c_complex(lua_State *L)                                        /* (-2,+1,e) */
{
    double x_re = luaL_checknumber(L, 1);
    double x_im = luaL_checknumber(L, 2);
    QLA_D_Complex *z = qlua_newComplex(L);
    
    QLA_real(*z) = x_re;
    QLA_imag(*z) = x_im;

    return 1;
}

static int
c_fmt(lua_State *L)                                            /* (-1,+1,e) */
{
    char c[72];
    QLA_D_Complex *z = qlua_checkComplex(L, 1);

    sprintf(c, "complex(%g, %g)", QLA_real(*z), QLA_imag(*z));
    lua_pushstring(L, c);

    return 1;
}

static int
c_re(lua_State *L)                                            /* (-1,+1,e) */
{
    QLA_D_Complex *z = qlua_checkComplex(L, 1);

    lua_pushnumber(L, QLA_real(*z));

    return 1;
}

static int
c_im(lua_State *L)                                            /* (-1,+1,e) */
{
    QLA_D_Complex *z = qlua_checkComplex(L, 1);
    
    lua_pushnumber(L, QLA_imag(*z));

    return 1;
}

static int
c_conj(lua_State *L)                                            /* (-1,+1,e) */
{
    QLA_D_Complex *z = qlua_checkComplex(L, 1);
    QLA_D_Complex *q = qlua_newComplex(L);

    QLA_real(*q) = QLA_real(*z);
    QLA_imag(*q) = -QLA_imag(*z);

    return 1;
}

static int
c_abs(lua_State *L)                                            /* (-1,+1,e) */
{
    QLA_D_Complex *z = qlua_checkComplex(L, 1);

    lua_pushnumber(L, hypot(QLA_real(*z), QLA_imag(*z)));

    return 1;
}

static int
c_neg(lua_State *L)                                            /* (-1,+1,e) */
{
    QLA_D_Complex *z = qlua_checkComplex(L, 1);
    QLA_D_Complex *q = qlua_newComplex(L);

    QLA_c_eqm_c(*q, *z);

    return 1;
}

static int
q_r_add_c(lua_State *L)                             /* (-2,+1,-) */
{
    QLA_D_Real a = luaL_checknumber(L, 1);
    QLA_D_Complex *b = qlua_checkComplex(L, 2);
    QLA_D_Complex *c = qlua_newComplex(L);

    QLA_c_eq_c(*c, *b);
    QLA_real(*c) += a;

    return 1;
}

static int
q_c_add_r(lua_State *L)                             /* (-2,+1,-) */
{
    QLA_D_Complex *b = qlua_checkComplex(L, 1);
    QLA_D_Real a = luaL_checknumber(L, 2);
    QLA_D_Complex *c = qlua_newComplex(L);

    QLA_c_eq_c(*c, *b);
    QLA_real(*c) += a;

    return 1;
}

static int
q_c_add_c(lua_State *L)                             /* (-2,+1,-) */
{
    QLA_D_Complex *a = qlua_checkComplex(L, 1);
    QLA_D_Complex *b = qlua_checkComplex(L, 2);
    QLA_D_Complex *c = qlua_newComplex(L);

    QLA_c_eq_c_plus_c(*c, *a, *b);

    return 1;
}

static int
q_r_sub_c(lua_State *L)                             /* (-2,+1,-) */
{
    QLA_D_Real a = luaL_checknumber(L, 1);
    QLA_D_Complex *b = qlua_checkComplex(L, 2);
    QLA_D_Complex *c = qlua_newComplex(L);

    QLA_c_eqm_c(*c, *b);
    QLA_real(*c) += a;

    return 1;
}

static int
q_c_sub_r(lua_State *L)                             /* (-2,+1,-) */
{
    QLA_D_Complex *b = qlua_checkComplex(L, 1);
    QLA_D_Real a = luaL_checknumber(L, 2);
    QLA_D_Complex *c = qlua_newComplex(L);

    QLA_c_eq_c(*c, *b);
    QLA_real(*c) -= a;

    return 1;
}

static int
q_c_sub_c(lua_State *L)                             /* (-2,+1,-) */
{
    QLA_D_Complex *a = qlua_checkComplex(L, 1);
    QLA_D_Complex *b = qlua_checkComplex(L, 2);
    QLA_D_Complex *c = qlua_newComplex(L);

    QLA_c_eq_c_minus_c(*c, *a, *b);

    return 1;
}

static int
q_r_mul_c(lua_State *L)                             /* (-2,+1,-) */
{
    QLA_D_Real a = luaL_checknumber(L, 1);
    QLA_D_Complex *b = qlua_checkComplex(L, 2);
    QLA_D_Complex *c = qlua_newComplex(L);

    QLA_c_eq_c_times_r(*c, *b, a);

    return 1;
}

static int
q_c_mul_r(lua_State *L)                             /* (-2,+1,-) */
{
    QLA_D_Complex *b = qlua_checkComplex(L, 1);
    QLA_D_Real a = luaL_checknumber(L, 2);
    QLA_D_Complex *c = qlua_newComplex(L);

    QLA_c_eq_c_times_r(*c, *b, a);

    return 1;
}

static int
q_c_mul_c(lua_State *L)                             /* (-2,+1,-) */
{
    QLA_D_Complex *a = qlua_checkComplex(L, 1);
    QLA_D_Complex *b = qlua_checkComplex(L, 2);
    QLA_D_Complex *c = qlua_newComplex(L);

    QLA_c_eq_c_times_c(*c, *a, *b);

    return 1;
}

static int
q_r_div_c(lua_State *L)                              /* (-2,+1,-) */
{
    double n;
    QLA_D_Real a = luaL_checknumber(L, 1);
    QLA_D_Complex *b = qlua_checkComplex(L, 2);
    QLA_D_Complex *c = qlua_newComplex(L);
    QLA_D_Complex bx;

    n = a / QLA_norm2_c(*b);
    QLA_c_eq_ca(bx, *b);
    
    QLA_c_eq_c_times_r(*c, bx, n);

    return 1;
}

static int
q_c_div_r(lua_State *L)                              /* (-2,+1,-) */
{
    QLA_D_Complex *b = qlua_checkComplex(L, 1);
    QLA_D_Real a = luaL_checknumber(L, 2);
    QLA_D_Complex *c = qlua_newComplex(L);

    QLA_c_eq_c_div_r(*c, *b, a);

    return 1;
}

static int
q_c_div_c(lua_State *L)                             /* (-2,+1,-) */
{
    QLA_D_Complex *a = qlua_checkComplex(L, 1);
    QLA_D_Complex *b = qlua_checkComplex(L, 2);
    QLA_D_Complex *c = qlua_newComplex(L);

    QLA_c_eq_c_div_c(*c, *a, *b);

    return 1;
}

int
q_c_gaussian(lua_State *L)
{
    mSeqRandom *s = qlua_checkSeqRandom(L, 1);
    QLA_D_Complex *z = qlua_newComplex(L);

    QLA_D_C_eq_gaussian_S(z, &s->state);

    return 1;
}

static const luaL_Reg mtComplex[] = {
    { "__tostring", c_fmt },
    { "__newindex", qlua_nowrite },
    { "real",       c_re },
    { "imag",       c_im },
    { "conj",       c_conj },
    { "abs",        c_abs },
    { "__unm",      c_neg },
    { "__add",      qlua_add },
    { "__sub",      qlua_sub },
    { "__mul",      qlua_mul },
    { "__div",      qlua_div },
    { NULL,         NULL }
};

static const luaL_Reg fComplex[] = {
    { "complex",   c_complex },
    { NULL,        NULL }
};

int
init_complex(lua_State *L)
{
    static const QLUA_Op2 ops[] = {
        {qlua_add_table, qReal,    qComplex, q_r_add_c},
        {qlua_add_table, qComplex, qReal,    q_c_add_r},
        {qlua_add_table, qComplex, qComplex, q_c_add_c},
        {qlua_sub_table, qReal,    qComplex, q_r_sub_c},
        {qlua_sub_table, qComplex, qReal,    q_c_sub_r},
        {qlua_sub_table, qComplex, qComplex, q_c_sub_c},
        {qlua_mul_table, qReal,    qComplex, q_r_mul_c},
        {qlua_mul_table, qComplex, qReal,    q_c_mul_r},
        {qlua_mul_table, qComplex, qComplex, q_c_mul_c},
        {qlua_div_table, qReal,    qComplex, q_r_div_c},
        {qlua_div_table, qComplex, qReal,    q_c_div_r},
        {qlua_div_table, qComplex, qComplex, q_c_div_c},
        {NULL,           qOther,   qOther,   NULL}
    };

    lua_getglobal(L, "_G");
    luaL_register(L, NULL, fComplex);
    qlua_metatable(L, mtnComplex, mtComplex, qComplex);
    qlua_reg_op2(ops);

    return 0;
}

int
fini_complex(lua_State *L)
{
    return 0;
}
