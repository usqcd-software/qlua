#ifndef MARK_F1F13269_3CC9_4F26_8140_EF8422F7ADE6
#define MARK_F1F13269_3CC9_4F26_8140_EF8422F7ADE6

enum {
    qG_z,  /* zero */
    qG_p,  /* +1 */
    qG_m,  /* -1 */
    qG_r,  /* real */
    qG_c,  /* complex */
    qG_t
};

typedef struct mGamma_s {
    int t;
    QLA_D_Real r;
    QLA_D_Complex c;
} mGamma;

typedef struct mClifford_s {
    mGamma g[16];
} mClifford;

mClifford *qlua_checkClifford(lua_State *L, int idx);
#ifdef HAS_GSL
#include "qmatrix.h"
mMatComplex *gamma2matrix(lua_State *L, int idx);
#endif

int init_gamma(lua_State *L);
void fini_gamma(void);

#endif /* !defined(MARK_F1F13269_3CC9_4F26_8140_EF8422F7ADE6) */
