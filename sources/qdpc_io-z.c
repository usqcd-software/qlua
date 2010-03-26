/* Readers and writers for a given colored type
 * expects
 * #define Qs(a)
 * #define Qx(a,b)
 * #define QC(x)
 * #define Qtype
 * #define Qop(x)
 * #define QNc
 * #define Qcolors
 *
 * Used with the following values:
 *    Qs      Qx            QC       Qtype       Qop             QNc   Qcolors
 *   QN(a,2)  QN(a,2) ## b  2        QT(QDP_D2)  QDP_D2_##x##_   '2'   "2"
 *   QN(a,3)  QN(a,2) ## b  3        QT(QDP_D3)  QDP_D3_##x##_   '3'   "3"
 *   QN(a,N)  QN(a,N) ## b  (x)->nc  QT(QDP_DN)  QDP_DN_##x##_   'N'   "N"
 */
static int
Qs(r_)(lua_State *L, mLattice *S, int Sidx, mReader *reader, int off, int nc)
{
    check_reader(L, reader);

    /* collect garbage now */
    CALL_QDP(L);
    /* read either one or n lattice color matrices */
    switch (lua_gettop(L) - off) {
    case 2: {
        /* one colored object */
        QDP_String *info = QDP_string_create();
        Qs(m) *X = Qs(qlua_new)(L, Sidx, nc);
        int status = Qop(read)(reader->ptr, info, X->ptr);
        if (status == 0) {
            /* read successful -- convert to LUA */
            lua_pushstring(L, QDP_string_ptr(info));
            QDP_string_destroy(info);
            return 2;
        }
        /* read failed */
        QDP_string_destroy(info);

        break;
    }
    case 3: {
        /* array of argv[1] color matrices */
        int n = luaL_checkint(L, 2 + off);
        int i;
        QDP_String *info;
        int status;
        
        /* sanity check */
        if (n <= 0)
            return luaL_error(L, "field count out of range");

        Qtype *X[n];
        lua_createtable(L, n, 0);
        for (i = 0; i < n; i++) {
            Qs(m) *xi = Qs(qlua_new)(L, Sidx, nc);
            X[i] = xi->ptr;
            lua_rawseti(L, -2, i + 1);
        }
        info = QDP_string_create();

        /* do the reading */
        status = Qop(vread)(reader->ptr, info, X, n);
        if (status == 0) {
            lua_pushstring(L, QDP_string_ptr(info));
            QDP_string_destroy(info);
            return 2;
        }
        /* read failed */
        QDP_string_destroy(info);
        for (i = 0; i < n; i++)
            Qop(destroy)(X[i]);

        break;
    }
    }
    return luaL_error(L, "qdpc read error");
}

static int
Qs(w_)(lua_State *L, mLattice *S, int Sidx, mWriter *writer)
{
    Qs(m) *X = Qs(qlua_check)(L, 2, S, -1);
    const char *info = luaL_checkstring(L, 3);
    QDP_String *xml;
    int status;

    check_writer(L, writer);

    /* collect garbage */
    CALL_QDP(L);

    /* prepare info string */
    xml = QDP_string_create();
    QDP_string_set(xml, (char *)info); /* [ sic ] */

    status = Qop(write)(writer->ptr, xml, X->ptr);
    QDP_string_destroy(xml);
    if (status == 0) {
        /* success -- return true */
        lua_pushboolean(L, 1);

        return 1;
    }
    return luaL_error(L, "qdpc write error");
}

static int
Qs(wt_)(lua_State *L, mLattice *S, int Sidx, mWriter *writer, int nc)
{
    int n = lua_objlen(L, 2);
    const char *info = luaL_checkstring(L, 3);
    int i;
    int status;
    QDP_String *xml;

    check_writer(L, writer);

    if (n <= 0)
        return luaL_error(L, "qdpc.write: bad table ");

    /* collect garbage */
    CALL_QDP(L);

    /* prepare info string */
    xml = QDP_string_create();
    QDP_string_set(xml, (char *)info); /* [ sic ] */

    Qtype *X[n];
    for (i = 0; i < n; i++) {
        Qs(m) *xi;
        /* full table indexing here */
        lua_pushnumber(L, i + 1); /* [ sic ] lua indexing */
        lua_gettable(L, 2);
        xi = Qs(qlua_check)(L, -1, S, nc);
        X[i] = xi->ptr;
        }
        
    /* do the write */
    status = Qop(vwrite)(writer->ptr, xml, X, n);
    QDP_string_destroy(xml);
    if (status == 0) {
        /* success -- clean up everything and return true */
        lua_pushboolean(L, 1);

        return 1;
    }
    return luaL_error(L, "qdpc write error");
}

#undef Qs
#undef Qx
#undef QC
#undef QNc
#undef Qop
#undef Qtype
#undef Qcolors
