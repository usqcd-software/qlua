/* Readers and writers for a given type (for all values of color)
 * expects
 *
 * #define QT(a)
 * #define QTx(a,b)
 * #define QN(a,b)
 * #define QA(p)
 *
 * Used for the following Lattice objects:
 *  QT                    QTx                      QN                  QA
 *  a##ColorVector        a##ColorVector##b        a##LatColVec##b     p##V
 *  a##ColorMatrix        a##ColorMatrix##b        a##LatColMat##b     p##M
 *  a##DiracFermion       a##DiracFermion##b       a##LatDirFerm##b    p##D
 *  a##DiracPropagator    a##DiracPropagator##b    a##LatDirProp##b    p##P
 */

#if USE_Nc2
#define Qs(a)    QN(a,2)
#define Qx(a,b)  QN(a,2) ## b
#define QC(x)    2
#define Qtype    QT(QDP_D2_)
#define Qop(x)   QA(QDP_D2_ ## x ## _)
#define QNc     '2'
#define Qcolors  "2"
#include "qdpc_io-z.c"                                               /* DEPS */
#endif

#if USE_Nc3
#define Qs(a)    QN(a,3)
#define Qx(a,b)  QN(a,3) ## b
#define QC(x)    3
#define Qtype    QT(QDP_D3_)
#define Qop(x)   QA(QDP_D3_ ## x ## _)
#define QNc     '3'
#define Qcolors  "3"
#include "qdpc_io-z.c"                                               /* DEPS */
#endif

#if USE_NcN
#define Qs(a)    QN(a,N)
#define Qx(a,b)  QN(a,N) ## b
#define QC(x)    (x)->nc
#define Qtype    QT(QDP_DN_)
#define Qop(x)   QA(QDP_DN_ ## x ## _)
#define QNc     'N'
#define Qcolors  "N"
#include "qdpc_io-z.c"                                               /* DEPS */
#endif

static int
QT(qdpc_r_)(lua_State *L)
{
#if USE_Nc2 || USE_Nc3 || USE_NcN
    mReader *reader = q_checkReader(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
#else
    mLattice *S = qlua_ObjLattice(L, 1);
#endif

    switch (S->nc) {
#if USE_Nc2
    case 2:  return QN(r_,2)(L, S, Sidx, reader, 0, 2);
#endif
#if USE_Nc3
    case 3:  return QN(r_,3)(L, S, Sidx, reader, 0, 3);
#endif
#if USE_NcN
    default: return QN(r_,N)(L, S, Sidx, reader, 0, S->nc);
#else
    default: return luaL_error(L, "bad number of colors");
#endif
    }
}

static int
QTx(qdpc_r_,N)(lua_State *L)
{
#if USE_Nc2 || USE_Nc3 || USE_NcN
    mReader *reader = q_checkReader(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);
#endif
    int nc = luaL_checkint(L, 2);

    switch (nc) {
#if USE_Nc2
    case 2:  return QN(r_,2)(L, S, Sidx, reader, 1, 2);
#endif
#if USE_Nc3
    case 3:  return QN(r_,3)(L, S, Sidx, reader, 1, 3);
#endif
#if USE_NcN
    default: return QN(r_,N)(L, S, Sidx, reader, 1, nc);
#else
    default: return luaL_error(L, "bad number of colors");
#endif
    }
}

static int
QT(qdpc_w_)(lua_State *L)
{
    mWriter *writer = q_checkWriter(L, 1, NULL);
    mLattice *S = qlua_ObjLattice(L, 1);
    int Sidx = lua_gettop(L);

    switch (qlua_qtype(L, 2)) {
#if USE_Nc2
    case QN(q,2): return QN(w_,2)(L, S, Sidx, writer);
#endif
#if USE_Nc3
    case QN(q,3): return QN(w_,3)(L, S, Sidx, writer);
#endif
#if USE_NcN
    case QN(q,N): return QN(w_,N)(L, S, Sidx, writer);
#endif
    case qTable: {
        lua_pushnumber(L, 1);
        lua_gettable(L, 2);
        switch (qlua_qtype(L, -1)) {
#if USE_Nc2
        case QN(q,2): return QN(wt_,2)(L, S, Sidx, writer, 2);
#endif
#if USE_Nc3
        case QN(q,3): return QN(wt_,3)(L, S, Sidx, writer, 3);
#endif
#if USE_NcN
        case QN(q,N): {
            QN(m,N) *X = QN(qlua_check,N)(L, -1, S, -1);
            return QN(wt_,N)(L, S, Sidx, writer, X->nc);
        }
#endif
        default:
            break;
        }
        break;
    }
    default:
        break;
    }
    return luaL_error(L, "qdpc write error");
}

#undef QT
#undef QTx
#undef QN
#undef QA
