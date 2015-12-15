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
