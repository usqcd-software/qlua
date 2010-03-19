#ifndef MARK_C8FCBCFB_8280_4F1D_9159_8E1C3F95E045
#define MARK_C8FCBCFB_8280_4F1D_9159_8E1C3F95E045

typedef enum {
    qss_all,
    qss_even,
    qss_odd,
    qss_dynamic,
    qss_slice,
    qss_upper,
    qss_lower
} qSubsetClass;

typedef struct {
    qSubsetClass  cl;
    int           axis;
    int           position;
} mLatSubset;

int init_latsubset(lua_State *L);
int fini_latsubset(lua_State *L);

#endif /* !defined(MARK_C8FCBCFB_8280_4F1D_9159_8E1C3F95E045) */
