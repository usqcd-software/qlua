/* Readers and writers for a given type
 * expects
 *  #define T_QTYPE  <name of the QDP type>
 *  #define QLUA_NAME(x)   x ## <QLUA type>
 *  #define X_ID(x)        x ## <QDP suffix>
 *  #define X_DF   / or #undef X_DF
 *  #define X_ID2(a,b)     a ## <QDP suffix> ## b ## <QDP suffix>
 *
 * Used for the following Lattice objects:
 *   T_QTYPE            QLUA_NAME         X_ID     X_DF
 *   QDP_Int            x ## LatInt       x ## I   undef
 *   QDP_RandomState    x ## LatRandom    x ## S   undef
 *   QDP_D_Real         x ## LatReal      x ## R   defined
 *   QDP_D_Complex      x ## LatComplex   x ## C   defined
 */

static int
X_ID(qdpc_r_)(lua_State *L)
{
    mReader *reader = q_checkReader(L, 1, NULL);
    int Sidx;

    check_reader(L, reader);
#ifdef X_DF
    mLattice *S = qlua_ObjLattice(L, 1);
#else
    qlua_ObjLattice(L, 1);
#endif
    Sidx = lua_gettop(L);

    /* collect garbage now */
    CALL_QDP(L);
    /* read either one or n lattice color matrices */
    switch (lua_gettop(L) - 1) {
    case 1: {
        /* one color matrix */
        QDP_String *info = QDP_string_create();
        QLUA_NAME(m) *U = QLUA_NAME(qlua_newZero)(L, Sidx);
        int status;
        
#ifdef X_DF
        if (get_prec(reader->ptr) == 'F') {
          T_sTYPE *tt = SopL(create)(S->lat);
          status = Sop(read)(reader->ptr, info, tt);
          FtoD(U->ptr, tt, S->all);
          Sop(destroy)(tt);
        } else {
#endif /* defined (X_DF) */
          status = X_ID(QDP_read_)(reader->ptr, info, U->ptr);
#ifdef X_DF
        }
#endif /* defined (X_DF) */

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
    case 2: {
        /* array of argv[1] color matrices */
        int n = luaL_checkint(L, 2);
        int i;
        QDP_String *info;
        int status;
        
        /* sanity check */
        if (n <= 0)
            return luaL_error(L, "field count out of range");

        T_QTYPE **U = qlua_malloc(L, n * sizeof (T_QTYPE *));
        lua_createtable(L, n, 0);
        for (i = 0; i < n; i++) {
            QLUA_NAME(m) *ui = QLUA_NAME(qlua_newZero)(L, Sidx);
            U[i] = ui->ptr;
            lua_rawseti(L, -2, i + 1);
        }
        info = QDP_string_create();

#ifdef X_DF
        if (get_prec(reader->ptr) == 'F') {
          T_sTYPE *tt[n];
          for (i = 0; i < n; i++) {
            tt[i] = SopL(create)(S->lat);
          }
          status = Sop(vread)(reader->ptr, info, tt, n);
          for (i = 0; i < n; i++) {
            FtoD(U[i], tt[i], S->all);
            Sop(destroy)(tt[i]);
          }
        } else {
#endif /* defined (X_DF) */
          /* do the reading */
          status = X_ID(QDP_vread_)(reader->ptr, info, U, n);
#ifdef X_DF
        }
#endif /* defined (X_DF) */
        if (status == 0) {
          lua_pushstring(L, QDP_string_ptr(info));
          QDP_string_destroy(info);
          qlua_free(L, U);
          return 2;
        }
        /* read failed */
        QDP_string_destroy(info);
        for (i = 0; i < n; i++)
            X_ID(QDP_destroy_)(U[i]);
        qlua_free(L, U);
        break;
    }
    }
    return luaL_error(L, "qdpc read error");
}

static int
X_ID(qdpc_w_)(lua_State *L)
{
    mWriter *writer = q_checkWriter(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    int status;
#ifdef X_DF
    int format = 'D';
#endif
    int Didx = 2;

    check_writer(L, writer);

#ifdef X_DF
    if (qlua_qtype(L, 2) == qString) {
        format = qlua_qio_file_precision(L, 2);
        Didx = 3;
    }
#endif
    const char *info = luaL_checkstring(L, Didx + 1);
    QDP_String *xml;

    /* collect garbage */
    CALL_QDP(L);

    /* prepare info string */
    xml = QDP_string_create();
    QDP_string_set(xml, (char *)info); /* [ sic ] */

    /* dispatch between single value and vector write */
    switch (qlua_qtype(L, Didx)) {
    case QLUA_NAME(q): {
        QLUA_NAME(m) *X = QLUA_NAME(qlua_check)(L, Didx, S);

#ifdef X_DF
        switch (format) {
        case 'F': {
            T_sTYPE *Y = X_ID2(QDP_F_create_,_L)(S->lat);
            if (Y == 0)
                return luaL_error(L, "not enough memory");
            X_ID3(QDP_FD_,_eq_)(Y, X->ptr, S->all);
            status = X_ID(QDP_F_write_)(writer->ptr, xml, Y);
            X_ID(QDP_F_destroy_)(Y);
            break;
        }
        case 'D':
            status = X_ID(QDP_D_write_)(writer->ptr, xml, X->ptr);
            break;
        default:
            return luaL_error(L, "unsupported QIO format");
        }
#else
        status = X_ID(QDP_write_)(writer->ptr, xml, X->ptr);
#endif
        if (status == 0) {
            /* success -- return true */
            QDP_string_destroy(xml);
            lua_pushboolean(L, 1);

            return 1;
        }
        break;
    }
    case qTable: {
        QLUA_NAME(m) *xi;
        int n, i;

        n = lua_objlen(L, Didx);
        if (n <= 0) {
            QDP_string_destroy(xml);
            return luaL_error(L, "qdpc.write: bad table ");
        }
        T_QTYPE **X = qlua_malloc(L, n * sizeof (T_QTYPE *));
        for (i = 0; i < n; i++) {
            /* full table indexing here */
            lua_pushnumber(L, i + 1); /* [ sic ] lua indexing */
            lua_gettable(L, Didx);
            xi = QLUA_NAME(qlua_check)(L, -1, S);
            X[i] = xi->ptr;
	    lua_pop(L, 1);
        }

#ifdef X_DF
        switch (format) {
        case 'F': {
            T_sTYPE **Y = qlua_malloc(L, n * sizeof (T_sTYPE *));
            for (i = 0; i < n; i++) {
                Y[i] = X_ID2(QDP_F_create_, _L)(S->lat);
                if (Y[i] == 0) {
                                        qlua_free(L, X);
                                        qlua_free(L, Y);
                    return luaL_error(L, "not enough memory");
                                }
                X_ID3(QDP_FD_,_eq_)(Y[i], X[i], S->all);
            }
            status = X_ID(QDP_F_vwrite_)(writer->ptr, xml, Y, n);
            for (i = 0; i < n; i++)
                X_ID(QDP_F_destroy_)(Y[i]);
                        qlua_free(L, Y);
            break;
        }
        case 'D':
            status = X_ID(QDP_D_vwrite_)(writer->ptr, xml, X, n);
            break;
        default:
                        qlua_free(L, X);
            return luaL_error(L, "unsupported QIO format");
        }
#else
        status = X_ID(QDP_vwrite_)(writer->ptr, xml, X, n);
#endif        
        if (status == 0) {
            /* success -- clean up everything and return true */
            QDP_string_destroy(xml);
            lua_pushboolean(L, 1);
                        qlua_free(L, X);

            return 1;
        }
        /* failure to write */
                qlua_free(L, X);
        break;
    }
    default:
        break;
    }

    /* cleanup */
    QDP_string_destroy(xml);

    return luaL_error(L, "qdpc write error");
}

#ifdef X_DF
#undef Sop
#undef SopL
#undef FtoD
#endif /* defined (X_DF) */
#undef T_QTYPE
#undef QLUA_NAME
#undef X_ID
#undef X_ID2
#undef X_ID3
#undef X_DF
#undef T_sTYPE
