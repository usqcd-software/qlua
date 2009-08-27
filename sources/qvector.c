#include <qlua.h>
#include <qcomplex.h>
#include <qvector.h>

const char *mtnVecInt     = "qcd.mtVecInt";
const char *mtnVecDouble  = "qcd.mtVecDouble";
const char *mtnVecComplex = "qcd.mtVecComplex";

#define VSIZE(s,n,t) (sizeof (s) + ((n)-1)*sizeof (t))

static void *
q_checkVector(lua_State *L, int idx, const char *mt)
{
    void *v = luaL_checkudata(L, idx, mt);

    luaL_argcheck(L, v != NULL, idx, "qcd.vector expected");

    return v;
}

tVecInt *
q_checkVecInt(lua_State *L, int idx)
{
    return (tVecInt *)q_checkVector(L, idx, mtnVecInt);
}

tVecDouble *
q_checkVecDouble(lua_State *L, int idx)
{
    return (tVecDouble *)q_checkVector(L, idx, mtnVecDouble);
}

tVecComplex *
q_checkVecComplex(lua_State *L, int idx)
{
    return (tVecComplex *)q_checkVector(L, idx, mtnVecComplex);
}

static void *
q_newVector(lua_State *L, int s, const char *mt)
{
    void *v = lua_newuserdata(L, s);

    luaL_getmetatable(L, mt);
    lua_setmetatable(L, -2);
    
    return v;
}

tVecInt *
q_newVecInt(lua_State *L, int size)
{
    return q_newVector(L, VSIZE(tVecInt, size, QLA_Int),
                       mtnVecInt);
}

tVecDouble *
q_newVecDouble(lua_State *L, int size)
{
    return q_newVector(L, VSIZE(tVecDouble, size, QLA_D_Real),
                       mtnVecDouble);
}

tVecComplex *
q_newVecComplex(lua_State *L, int size)
{
    return q_newVector(L, VSIZE(tVecComplex, size, QLA_D_Complex),
                       mtnVecComplex);
}

/* integer vectors */
static int
vi_fmt(lua_State *L)                                           /* (-1,+1,e) */
{
    char fmt[72];
    tVecInt *v = q_checkVecInt(L, 1);

    sprintf(fmt, "vec_int[%d]", v->size);
    lua_pushstring(L, fmt);

    return 1;
}

static int
vi_get(lua_State *L)                                           /* (-2,+1,e) */
{
    tVecInt *v = q_checkVecInt(L, 1);
    
    switch(qlua_gettype(L, 2)) {
    case qString: {
        lua_pushstring(L, "length");
        if (lua_equal(L, 2, -1)) {
            lua_pushnumber(L, v->size);
            return 1;
        }
    }
    case qReal: {
        int i = luaL_checkint(L, 2);
        if ((i >= 0) && (i < v->size)) {
            lua_pushnumber(L, v->val[i]);
            return 1;
        }
    }
    }
    return luaL_error(L, "bad index for vec_int[]");
}

static int
vi_put(lua_State *L)                                           /* (-3,+0,e) */
{
    tVecInt *v = q_checkVecInt(L, 1);
    int idx = luaL_checkint(L, 2);
    int x = luaL_checkint(L, 3);

    if ((idx >= 0) && (idx < v->size)) {
        v->val[idx] = x;
        return 0;
    }
    return luaL_error(L, "bad index for vec_int[]");
}

static int
v_int(lua_State *L)                                            /* (-1,+1,e) */
{
    int s = luaL_checkint(L, 1);
    tVecInt *v = q_newVecInt(L, s);
    v->size = s;
    
    return 1;
}

static const luaL_Reg mtVecInt[] = {
    { "__tostring",     vi_fmt    },
    { "__index",        vi_get    },
    { "__newindex",     vi_put    },
    { NULL,             NULL      }
};

/* double vectors */
static int
vd_fmt(lua_State *L)                                           /* (-1,+1,e) */
{
    char fmt[72];
    tVecDouble *v = q_checkVecDouble(L, 1);

    sprintf(fmt, "vec_double[%d]", v->size);
    lua_pushstring(L, fmt);

    return 1;
}

static int
vd_get(lua_State *L)                                           /* (-2,+1,e) */
{
    tVecDouble *v = q_checkVecDouble(L, 1);
    
    switch(qlua_gettype(L, 2)) {
    case qString: {
        lua_pushstring(L, "length");
        if (lua_equal(L, 2, -1)) {
            lua_pushnumber(L, v->size);
            return 1;
        }
    }
    case qReal: {
        int i = luaL_checkint(L, 2);
        if ((i >= 0) && (i < v->size)) {
            lua_pushnumber(L, v->val[i]);
            return 1;
        }
    }
    }
    return luaL_error(L, "bad index for vec_double[]");
}

static int
vd_put(lua_State *L)                                           /* (-3,+0,e) */
{
    tVecDouble *v = q_checkVecDouble(L, 1);
    int idx = luaL_checkint(L, 2);
    double x = luaL_checknumber(L, 3);

    if ((idx >= 0) && (idx < v->size)) {
        v->val[idx] = x;
        return 0;
    }
    return luaL_error(L, "bad index for vec_double[]");
}

static int
v_double(lua_State *L)                                         /* (-1,+1,e) */
{
    int s = luaL_checkint(L, 1);
    tVecDouble *v = q_newVecDouble(L, s);
    v->size = s;
    
    return 1;
}

static const luaL_Reg mtVecDouble[] = {
    { "__tostring",     vd_fmt    },
    { "__index",        vd_get    },
    { "__newindex",     vd_put    },
    { NULL,             NULL      }
};

/* complex vectors */
static int
vc_fmt(lua_State *L)                                           /* (-1,+1,e) */
{
    char fmt[72];
    tVecComplex *v = q_checkVecComplex(L, 1);

    sprintf(fmt, "vec_complex[%d]", v->size);
    lua_pushstring(L, fmt);

    return 1;
}

static int
vc_get(lua_State *L)                                           /* (-2,+1,e) */
{
    tVecComplex *v = q_checkVecComplex(L, 1);
    
    switch(qlua_gettype(L, 2)) {
    case qString: {
        lua_pushstring(L, "length");
        if (lua_equal(L, 2, -1)) {
            lua_pushnumber(L, v->size);
            return 1;
        }
    }
    case qReal: {
        int i = luaL_checkint(L, 2);
        if ((i >= 0) && (i < v->size)) {
            QLA_D_Complex *z = q_newComplex(L);
            QLA_c_eq_c(*z, v->val[i]);
            return 1;
        }
    }
    }

    return luaL_error(L, "bad index for vec_complex[]");
}

static int
vc_put(lua_State *L)                                           /* (-3,+0,e) */
{
    tVecComplex *v = q_checkVecComplex(L, 1);
    int idx = luaL_checkint(L, 2);

    if ((idx >= 0) && (idx < v->size)) {
        switch (qlua_gettype(L, 3)) {
        case qReal: {
            double x = luaL_checknumber(L, 3);
            QLA_real(v->val[idx]) = x;
            QLA_imag(v->val[idx]) = 0;
            return 0;
        }
        case qComplex: {
            QLA_D_Complex *z = q_checkComplex(L, 3);
            QLA_c_eq_c(v->val[idx], *z);
            return 0;
        }
        }
    }

    return luaL_error(L, "bad index for vec_double[]");
}

static int
v_complex(lua_State *L)                                         /* (-1,+1,e) */
{
    int s = luaL_checkint(L, 1);
    tVecComplex *v = q_newVecComplex(L, s);
    v->size = s;
    
    return 1;
}

static const luaL_Reg mtVecComplex[] = {
    { "__tostring",     vc_fmt    },
    { "__index",        vc_get    },
    { "__newindex",     vc_put    },
    { NULL,             NULL      }
};

/* vector constructors */
static const luaL_Reg fVector[] = {
    { "vec_int",       v_int     },
    { "vec_double",    v_double  },
    { "vec_complex",   v_complex },
    { NULL,            NULL      }
};

int
init_vector(lua_State *L)
{
    luaL_register(L, qcdlib, fVector);
    qlua_metatable(L, mtnVecInt,     mtVecInt);
    qlua_metatable(L, mtnVecDouble,  mtVecDouble);
    qlua_metatable(L, mtnVecComplex, mtVecComplex);
    return 0;
}

int
fini_vector(lua_State *L)
{
    return 0;
}
