#ifndef QLOPT_H_
#define QLOPT_H_

#include <complex.h>
#include <stddef.h>
#include <limits.h>
#include <float.h>

#include "qlua.h"

/* constants for qloptElem.key */
#define QLOPT_KEY  (-10000)        /* use key instead of position */
#define QLOPT_NEXT (-10001)        /* position next relative to the prev item; =1 if the first item in the table */
#define QLOPT_END  (0)             /* last element */

/* cases, in order of decreasing importance :
   1) mandatory with correct type : exit if missing or mismatch, print error + parsing stack; >> print warning
   2) optional, with correct type : exit, print error + parsing stack >> print warning
   3) missing warning 
   */

#define QLOPT_OPT_WARN     0x1
#define QLOPT_OPT          0x2 
#define QLOPT_BAD_WARN     0x4 /* if set, print warning, skip (cannot print stack) */
#define QLOPT_BAD_SKIP     0x8 /* if not set, exit + print error */

typedef int qlopt_Int;
typedef double qlopt_Real;
typedef complex double qlopt_Complex;



typedef enum qloptType_e {
    QLOPT_NONE = 0,     /* in parsing : skip ; in writing : put nil */
    /* scalar types */
    QLOPT_INT  = 1,
    QLOPT_BOOL,
    QLOPT_REAL,
    QLOPT_COMPLEX,
    QLOPT_STRING,

    /* lattice fields : copy ptr only, not the data
       (and hope it remains on stack long enough) */
    QLOPT_LATINT,
    QLOPT_LATREAL,
    QLOPT_LATCOMPLEX,
    QLOPT_LATCOLVEC,
    QLOPT_LATCOLMAT,
    QLOPT_LATDIRFERM,
    QLOPT_LATDIRPROP,
    
    /* other Qlua objects */
    QLOPT_GAMMA,
    QLOPT_SUBSET,

    /* nested collections */
    QLOPT_TABLE,       /* general table ; specify parsing by position/key */
    QLOPT_ARRAY,       /* uniform multidim.array of qloptType 
                           except QLOPT_TABLE, QLOPT_ARRAY, QLOPT_CALLBACK;
                           gnore position/key ; allocate if necessary */
    /* general */
//    QLOPT_CALLBACK_e
} qloptType;

/* overset of QLUA_Type suitable for controlling parsing options */
//typedef struct qloptType_s {
//    QLUA_Type       qt;
//    qloptOthType_e  qlopt_t;
//} qloptType;
//
//static qloptType QLOPT_NONE = { qOther, QLOPT_NONE_e } ;
//static qloptType QLOPT_INT = { qOther, QLOPT_INT_e } ;
//static qloptType QLOPT_BOOL = { qOther, QLOPT_BOOL_e } ;
//static qloptType QLOPT_REAL = { qReal, } ;
//static qloptType QLOPT_COMPLEX = { qOther, } ;
//static qloptType QLOPT_STRING = { qOther, } ;
//
//static qloptType QLOPT_LATINT = { qOther, } ;
//static qloptType QLOPT_LATREAL = { qOther, } ;
//static qloptType QLOPT_LATCOMPLEX = { qOther, } ;
//static qloptType QLOPT_LATCOLVEC = { qOther, } ;
//static qloptType QLOPT_LATCOLMAT = { qOther, } ;
//static qloptType QLOPT_LATDIRFERM = { qOther, } ;
//static qloptType QLOPT_LATDIRPROP = { qOther, } ;
//static qloptType QLOPT_GAMMA = { qOther, } ;
//static qloptType QLOPT_SUBSET = { qOther, } ;
//
//static qloptType QLOPT_TABLE = { qOther, QLOPT_TABLE };
//static qloptType QLOPT_ARRAY = { qOther, QLOPT_ARRAY };
//

union  qloptData;
struct qloptNdarray;
struct qloptElem;
//struct qloptCallback ;


#if 0
typedef int (*qloptCallbackFunc_t)(int stack_pos, void *arg);
struct qloptCallback {
    qloptCallbackFunc_t func;
    void *arg;
};
typedef struct qloptCallback qloptCallback_s;
#endif

/* a union describing a source or destination for stack writing/parsing, resp. */
typedef union qloptData_u {
    /*  int, double, complex double : copy to (*(<type> *)p)
        string : *p <- luaL_checkstring
        QDP_FIELD : copy ptr to (QDP_FIELD *)(*p) 
        subset : ??
        table : parse table

        callback - not implemented
     */
    off_t           off;    /* data offset if base address is specified */
    void            *ptr;   /* data location */
    struct qloptNdarray_s   
                    *arr;   /* nd-array */  
    struct qloptElem_s      
                    *tab;   /* parse a table: list of qloptElem_s terminated by QLOPT_END */
//   qloptCallback_s f;    /* parse with a callback */
#if 1
    /* typed aliases to void *p */
    qlopt_Int       *i;     /* int or bool */
    qlopt_Real      *r;     /* floating-point */
    qlopt_Complex   *c;

    const char      **s;    /* string is not allocated; only a pointer is copied */
#endif  

} qloptData;

/* n-dimensional arrays of reg.types or tables 
   XXX ALL FIELDS MUST BE INITIALIZED, OR YOUR CODE WILL DIE SCREAMING
   perhaps TODO init functions for different cases that will check all parameters for correctness
 */
typedef struct qloptNdarray_s {
    /* parse and store a multidimensional array of dim[0.. (ndim-1)] dimensions, 
       the last dim changes fastest (C-ordering) */
    int         ndim;       /* number of dimensions */
    int         *dim;       /* verify if dimensions >0, parse if -1 and allocate buf */
    int         *maxdim;    /* check that dimensions do not exceed the limit (ignore if NULL or entries <=0) */
    qloptType   argtype;    /* cannot be callback, ndarray 
                               TODO develop interface for ndarray of tables */
    void        *buf;       /* allocate if buf == NULL; report error if have undef(==-1) dim[] */
    size_t      argsize;    /* size (bytes) to allocate for 1 element; used only if buf==NULL */
    void        *fill;      /* if not NULL, prefill with the value [fill:fill+argsize) */
    struct qloptElem_s *elem_list; /* meaningful only if argtype == QLOPT_TABLE */
} qloptNdarray;

typedef struct qloptElem_s {
    int         pos;        /* if QLOPT_KEY, use key */
    char        *key;       /* if NULL, default callback for undef keys */
    qloptType   argtype;    /* what to look for */
    int         flags;      /* how to handle optional / type mismatch / etc */
    union qloptData_u   val;
} qloptElem;
    
/* arg parsing functions */
int qlopt_read_array(lua_State *L, int pos, qloptNdarray *a);
int qlopt_read_table(lua_State *L, int pos, qloptElem *elem_list);
int qlopt_read_stack(lua_State *L, qloptElem *elem_list);
int q_test_qlopt(lua_State *L);
/* use these functions to initialize qloptNdarray; direct structure init is discouraged! */
qloptNdarray qlopt_array_scalar(int ndim, int *dim, int *maxdim, 
        qloptType argtype, void *buf, void *fill);
qloptNdarray qlopt_array_struct(int ndim, int *dim, int *maxdim, 
        qloptElem *elem_list, size_t argsize, void *buf, void *fill);

#endif/*QLOPT_H_*/
