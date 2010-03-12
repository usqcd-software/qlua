/* Readers and writers for a given type
 * expects
 *  #define T_QTYPE  <name of the QDP type>
 *  #define QLUA_NAME(x)   x ## <QLUA type>
 *  #define X_ID(x)        x ## <QDP suffix>
 */

static int
X_ID(qdpc_r_)(lua_State *L)
{
    mReader *reader = q_checkReader(L, 1);

    check_reader(L, reader);

    /* collect garbage now */
    CALL_QDP(L);
    /* read either one or n lattice color matrices */
    switch (lua_gettop(L)) {
    case 1: {
        /* one color matrix */
        QDP_String *info = QDP_string_create();
        QLUA_NAME(m) *U = QLUA_NAME(qlua_new)(L);
        
        if (X_ID(QDP_read_)(reader->ptr, info, U->ptr) == 0) {
            /* read successful -- convert to LUA */
            lua_pushstring(L, QDP_string_ptr(info));
            QDP_string_destroy(info);
            return 2;
        }
        /* read failed */
        QDP_string_destroy(info);

        break;
    }
    case 2: {
        /* array of argv[1] color matrices */
        int n = luaL_checkint(L, 2);
        int i;
        QDP_String *info;
        T_QTYPE **U;
        int status;
        
        /* sanity check */
        if (n <= 0)
            return luaL_error(L, "field count out of range");

        U = qlua_malloc(L, n * sizeof (T_QTYPE *));
        lua_createtable(L, n, 0);
        for (i = 0; i < n; i++) {
            QLUA_NAME(m) *ui = QLUA_NAME(qlua_new)(L);
            U[i] = ui->ptr;
            lua_rawseti(L, -2, i + 1);
        }
        info = QDP_string_create();

        /* do the reading */
        status = X_ID(QDP_vread_)(reader->ptr, info, U, n);
        qlua_free(L, U);
        if (status == 0) {
            lua_pushstring(L, QDP_string_ptr(info));
            QDP_string_destroy(info);
            return 2;
        }
        /* read failed */
        QDP_string_destroy(info);
        for (i = 0; i < n; i++)
            X_ID(QDP_destroy_)(U[i]);
        
        break;
    }
    }
    return luaL_error(L, "qdpc read error");
}

static int
X_ID(qdpc_w_)(lua_State *L)
{
    mWriter *writer = q_checkWriter(L, 1);
    const char *info = luaL_checkstring(L, 3);
    QDP_String *xml;

    check_writer(L, writer);

    /* collect garbage */
    CALL_QDP(L);

    /* prepare info string */
    xml = QDP_string_create();
    QDP_string_set(xml, (char *)info); /* [ sic ] */

    /* dispatch between single value and vector write */
    switch (qlua_gettype(L, 2)) {
    case QLUA_NAME(q): {
        QLUA_NAME(m) *U = QLUA_NAME(qlua_check)(L, 2);

        if (X_ID(QDP_write_)(writer->ptr, xml, U->ptr) == 0) {
            /* success -- return true */
            QDP_string_destroy(xml);
            lua_pushboolean(L, 1);

            return 1;
        }
        break;
    }
    case qTable: {
        QLUA_NAME(m) *ui;
        T_QTYPE **U;
        int n, i;
        int status;

        n = lua_objlen(L, 2);
        if (n <= 0) {
            QDP_string_destroy(xml);
            return luaL_error(L, "qdpc.write: bad table ");
        }
        U = qlua_malloc(L, n * sizeof (T_QTYPE *));
        for (i = 0; i < n; i++) {
            /* full table indexing here */
            lua_pushnumber(L, i + 1); /* [ sic ] lua indexing */
            lua_gettable(L, 2);
            ui = QLUA_NAME(qlua_check)(L, -1);
            U[i] = ui->ptr;
        }
        
        /* do the write */
        status = X_ID(QDP_vwrite_)(writer->ptr, xml, U, n);
        qlua_free(L, U);
        if (status == 0) {
            /* success -- clean up everything and return true */
            QDP_string_destroy(xml);
            lua_pushboolean(L, 1);

            return 1;
        }
        /* failure to write */
        break;
    }
    }

    /* cleanup */
    QDP_string_destroy(xml);

    return luaL_error(L, "qdpc write error");
}

#undef T_QTYPE
#undef QLUA_NAME
#undef X_ID
