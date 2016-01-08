#ifndef QOPT_H_
#define QOPT_H_

#include <complex.h>
#include <stddef.h>
#include <limits.h>
#include <float.h>

#include "qlua.h"

/* constants for qoptElem.key */
#define QOPT_KEY  (-10000)        /* use key instead of position */
#define QOPT_NEXT (-10001)        /* position next relative to the prev item; =1 if the first item in the table */
#define QOPT_END  (0)             /* last element */
/* flags for handling missing or mismatching arguments */
#define QOPT_OPT_WARN     0x1
#define QOPT_OPT          0x2 
#define QOPT_BAD_WARN     0x4 /* if set, print warning, skip (cannot print stack) */
#define QOPT_BAD_SKIP     0x8 /* if not set, exit + print error */

typedef int qopt_Int;
typedef double qopt_Real;
typedef complex double qopt_Complex;


#if 0
typedef enum qoptType_e {
    qOptNone = 0,     /* in parsing : skip ; in writing : put nil */
    /* scalar types */
    qOptInt  = 1,
    qOptBool,
    qReal,
    qComplex,
    qString,

    /* lattice fields : copy ptr only, not the data
       (and hope it remains on stack long enough) */
    qLatInt,
    qLatReal,
    qLatComplex,
    qLatColVec3,
    qLatColMat3,
    qLatDirFerm3,
    qLatDirProp3,
    
    /* other Qlua objects */
    qGamma,
    qLatSubset,

    /* nested collections */
    qOptStruct,       /* general table ; specify parsing by position/key */
    qOptArray,       /* uniform multidim.array of qoptType 
                           except qOptStruct, qOptArray, QOPT_CALLBACK;
                           gnore position/key ; allocate if necessary */
    /* general */
//    QOPT_CALLBACK_e
} qoptType;
#endif

union  qoptData;
struct qoptNdarray;
struct qoptElem;
//struct qoptCallback ;


#if 0
typedef int (*qoptCallbackFunc_t)(int stack_pos, void *arg);
struct qoptCallback {
    qoptCallbackFunc_t func;
    void *arg;
};
typedef struct qoptCallback qoptCallback_s;
#endif

/* a union describing a source or destination for stack writing/parsing, resp. */
typedef union qoptData_u {
    /*  int, double, complex double : copy to (*(<type> *)p)
        string : *p <- luaL_checkstring
        QDP_FIELD : copy ptr to (QDP_FIELD *)(*p) 
        subset : ??
        table : parse table

        callback - not implemented
     */
    off_t           off;    /* data offset if base address is specified */
    void            *ptr;   /* data location */
    struct qoptNdarray_s   
                    *arr;   /* nd-array */  
    struct qoptElem_s      
                    *tab;   /* parse a table: list of qoptElem_s terminated by QOPT_END */
//   qoptCallback_s f;    /* parse with a callback */
#if 1
    /* typed aliases to void *p */
    qopt_Int       *i;     /* int or bool */
    qopt_Real      *r;     /* floating-point */
    qopt_Complex   *c;

    const char      **s;    /* string is not allocated; only a pointer is copied */
#endif  

} qoptData;

/* n-dimensional arrays of reg.types or tables 
   XXX ALL FIELDS MUST BE INITIALIZED, OR YOUR CODE WILL DIE SCREAMING
   perhaps TODO init functions for different cases that will check all parameters for correctness
 */
typedef struct qoptNdarray_s {
    /* parse and store a multidimensional array of dim[0.. (ndim-1)] dimensions, 
       the last dim changes fastest (C-ordering) */
    int         ndim;       /* number of dimensions */
    int         *dim;       /* verify if dimensions >0, parse if -1 and allocate buf */
    int         *maxdim;    /* check that dimensions do not exceed the limit (ignore if NULL or entries <=0) */
    QLUA_Type   argtype;    /* cannot be callback, ndarray 
                               TODO develop interface for ndarray of tables */
    void        *buf;       /* allocate if buf == NULL; report error if have undef(==-1) dim[] */
    size_t      argsize;    /* size (bytes) to allocate for 1 element; used only if buf==NULL */
    void        *fill;      /* if not NULL, prefill with the value [fill:fill+argsize) */
    struct qoptElem_s *elem_list; /* meaningful only if argtype == qOptStruct */
} qoptNdarray;

typedef struct qoptElem_s {
    int         pos;        /* if QOPT_KEY, use key */
    char        *key;       /* if NULL, default callback for undef keys */
    QLUA_Type   argtype;    /* what to look for */
    int         flags;      /* how to handle optional / type mismatch / etc */
    union qoptData_u   val;
} qoptElem;
    
/* arg parsing functions */
int qopt_read_array(lua_State *L, int pos, qoptNdarray *a);
int qopt_read_table(lua_State *L, int pos, qoptElem *elem_list);
int qopt_read_stack(lua_State *L, qoptElem *elem_list);
int q_test_qopt(lua_State *L);
/* use these functions to initialize qoptNdarray; direct structure init is discouraged! */
qoptNdarray qopt_array_scalar(int ndim, int *dim, int *maxdim, 
        QLUA_Type argtype, void *buf, void *fill);
qoptNdarray qopt_array_struct(int ndim, int *dim, int *maxdim, 
        qoptElem *elem_list, size_t argsize, void *buf, void *fill);

#endif/*QOPT_H_*/
