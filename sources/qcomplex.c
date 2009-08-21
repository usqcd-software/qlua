#include <qlua.h>
#include <math.h>
#include <qcomplex.h>

static const char *mt_name = "mtComplex";

typedef struct {
    double re, im;
} q_Complex;

static q_Complex *
q_checkComplex(lua_State *L, int idx)
{
    void *ud = luaL_checkudata(L, idx, mt_name);

    luaL_argcheck(L, ud != NULL, idx, "qcd.complex expected");

    return (q_Complex *)ud;
}

static q_Complex *
q_newComplex(lua_State *L)
{
    q_Complex *z = (q_Complex *)lua_newuserdata(L, sizeof (q_Complex));
    
    luaL_getmetatable(L, mt_name);
    lua_setmetatable(L, -2);

    return z;
}


static int
c_complex(lua_State *L)                                        /* (-2,+1,e) */
{
    double x_re = luaL_checknumber(L, 1);
    double x_im = luaL_checknumber(L, 2);
    q_Complex *z = q_newComplex(L);
    
    z->re = x_re;
    z->im = x_im;

    return 1;
}

static int
c_fmt(lua_State *L)                                            /* (-1,+1,e) */
{
    char c[72];
    q_Complex *z = q_checkComplex(L, 1);

    sprintf(c, "complex(%g, %g)", z->re, z->im);
    lua_pushstring(L, c);

    return 1;
}

static int
c_re(lua_State *L)                                            /* (-1,+1,e) */
{
    q_Complex *z = q_checkComplex(L, 1);

    lua_pushnumber(L, z->re);

    return 1;
}

static int
c_im(lua_State *L)                                            /* (-1,+1,e) */
{
    q_Complex *z = q_checkComplex(L, 1);

    lua_pushnumber(L, z->im);

    return 1;
}

static int
c_conj(lua_State *L)                                            /* (-1,+1,e) */
{
    q_Complex *z = q_checkComplex(L, 1);
    q_Complex *q = q_newComplex(L);

    q->re = z->re;
    q->im = -z->im;

    return 1;
}

static int
c_abs(lua_State *L)                                            /* (-1,+1,e) */
{
    q_Complex *z = q_checkComplex(L, 1);

    lua_pushnumber(L, hypot(z->re, z->im));

    return 1;
}

static int
c_neg(lua_State *L)                                            /* (-1,+1,e) */
{
    q_Complex *z = q_checkComplex(L, 1);
    q_Complex *q = q_newComplex(L);

    q->re = -z->re;
    q->im = -z->im;
    return 1;
}

static const luaL_Reg mtComplex[] = {
    { "__tostring", c_fmt },
    { "re",         c_re },
    { "im",         c_im },
    { "conj",       c_conj },
    { "abs",        c_abs },
    { "__unm",      c_neg },
#if 0 /* XXX */
    { "__add",      c_add },
    { "__sub",      c_sub },
    { "__mul",      c_mul },
    { "__div",      c_div },
#endif
    { NULL, NULL }
};

static const luaL_Reg fComplex[] = {
    { "complex",   c_complex },
    { NULL, NULL }
};

int
init_complex(lua_State *L)
{
    luaL_register(L, qcdlib, fComplex);
    qlua_metatable(L, mt_name, mtComplex);
    return 0;
}
