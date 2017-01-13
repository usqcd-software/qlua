#include "qopt.h"


/* aux: calculate a product of integer array */
static int 
prod_int_(int n, int *a) 
{
    int res;
    for (res = 1; n--; res *= *(a++));
    return res;
}

#define printf0(...) do { if (QDP_this_node == qlua_master_node) printf(__VA_ARGS__); } while(0)
static void
print_arr_general(int ndim, int *dim, void *buf, size_t argsize, 
        void (*print_rec)(void *buf, const char *msg_head),  /* print a record in 1 line */ 
        const char *msg_head)
{
    char strbuf[1024];
    if (0 == ndim)
        print_rec(buf, msg_head);
    else if (0 < ndim) {
        int memb_stride = prod_int_(ndim - 1, dim + 1);
        for (int i = 0 ; i < *dim ; i++) {
            snprintf(strbuf, sizeof(strbuf), "%s[%d]", msg_head, i);
            print_arr_general(ndim - 1, dim + 1, buf + memb_stride * argsize * i, 
                    argsize, print_rec, strbuf);
        }
    }
    /* if ndim < 0, do nothing */
}



static void
print_dim(int n, const int dim[], const char *msg_head)
{
    printf0("%s@{", msg_head);
    for (int i = 0 ; i < n ; i++)
        printf0("%d ", dim[i]);
    printf0("}\n");
}

/* example for reading an array of records */
typedef struct recx_s {
    QDP_D_Complex   *cplx;
    char            *name;
    QLA_Real        bc_t;
} recx;
qoptElem qopt_tab_recx[] = {
    { 1, NULL, qString,    0, {.off = offsetof(recx, name) } },
    { 2, NULL, qLatComplex,0, {.off = offsetof(recx, cplx) } },
    { 3, NULL, qReal,      0, {.off = offsetof(recx, bc_t) } },
    { QOPT_END }
};

static void
print_recx_(void *a_, const char *msg_head) 
{
    recx *a = (recx *)a_;
    QLA_D_Real n2 = 0.;
    QDP_D_r_eq_norm2_C(&n2, a->cplx, QDP_all_L(QDP_D_get_lattice_C(a->cplx)));
    printf0("%s = {name='%s'  f=%p(|f|2=%e)  bc_t=%f}\n", 
            msg_head, a->name, 
            (void*)a->cplx, (double)n2,
            (double)(a->bc_t));
}
static void
print_str_(void *a_, const char *msg_head)
{
    printf0("%s = '%s'\n", msg_head, *(char **)a_);
}
static void 
print_qoptint_(void *a_, const char *msg_head) 
{
    printf0("%s = %lld\n", msg_head, (long long int)(*(qopt_Int *)a_));
}
static void 
print_qoptreal_(void *a_, const char *msg_head) 
{
        printf0("%s = %e\n", msg_head, (double)(*(qopt_Real *)a_));
}
static void 
print_qoptcomplex_(void *a_, const char *msg_head) 
{
    qopt_Complex *a = (qopt_Complex *)a_;
    printf0("%s = (%e, %e)\n", msg_head, creal(*a), cimag(*a));
}


/* built-in test */
int 
q_test_qopt(lua_State *L)
{
    int i_arr[3];
    double d_arr[5];
    const char *s_arr[7];
    void *xx = NULL;

    /*
    qcd.test_qopt(
      i1, 
      i2, 
      r3, 
      c4,
      recx_MxN,
      { 
        cplxMxNxP,
        str2xN,
        struct1 = { s1i1, s1i2, s1r3, s1i4},
        struct2 = { s2i1, s2r2, s2r3, s2r4_opt },
        int2x3x4= int2x3x4,
        intNx4  = intNx4,
        realNx2 = realNx2,
        str1    = str1, 
        str2    = str2,
        bool1   = true
        }
        )
    */
    qopt_Int i1 = -1,
             i2 = -1;
    qopt_Real r3 = -1.;
    qopt_Complex c4 = -1.;


    qopt_Int s1i1  = -1,
             s1i2  = -1,
             s1i4  = -1;
    qopt_Real s1r3 = -1.;
    qoptElem qopt_tab_struct1[] = {
        { 1, NULL, qOptInt, 0,  { .i = &s1i1 } },
        { 2, NULL, qOptInt, 0,  { .i = &s1i2 } },
        { 3, NULL, qReal,   0,  { .r = &s1r3 } },
        { 4, NULL, qOptInt, 0,  { .i = &s1i4 } },
        { QOPT_END }
    };
    
    qopt_Int s2i1 = -1;
    qopt_Real s2r2  = DBL_MAX,
              s2r3  = DBL_MAX,
              s2r4  = DBL_MAX;
    qoptElem qopt_tab_struct2[] = {
        { 1, NULL, qOptInt, 0,  { .i = &s2i1 } },
        { 2, NULL, qReal,   0,  { .r = &s2r2 } },
        { 3, NULL, qReal,   0,  { .r = &s2r3 } },
        { 4, NULL, qReal,   QOPT_OPT, { .r = &s2r4 } },
        { QOPT_END }
    };
    
    int recx_dim[]     = { -1, -1 };
    int recx_dim_max[] = { -1, -1 };
    qoptArray qopt_arr_recx = qopt_array_struct(2, recx_dim, recx_dim_max, 
                qopt_tab_recx, sizeof(recx), NULL, NULL);

    int cplx_dim[]     = { -1, -1, -1 },
        cplx_dim_max[] = { 10, 10, 10 };
    qoptArray qopt_arr_cplx = qopt_array_scalar(3, cplx_dim, cplx_dim_max, 
                qComplex, NULL, NULL);

    int str_dim[]     = {  2, -1 },
        str_dim_max[] = {  4,  4 };
    const char *str_arr_buf[16];
    qoptArray qopt_arr_str = qopt_array_scalar(2, str_dim, str_dim_max, 
                qString, str_arr_buf, NULL);
    
    int int1_dim[] = {-1, 4},
        int1_dim_max[] = { 10, -1};
    qoptArray qopt_arr_int1 = qopt_array_scalar(2, int1_dim, int1_dim_max, 
                qOptInt, NULL, NULL);

    int int2_dim[] = { 2, 3, 4 };
    qopt_Int int2[2][3][4];
    qoptArray qopt_arr_int2 = qopt_array_scalar(3, int2_dim, NULL, 
                qOptInt, int2, NULL);
    
    int real_dim[]    = { -1, 2 },
        real_dim_max[] = { 10, -1 }; 
    qoptArray qopt_arr_real = { 2, real_dim, real_dim_max, 
                qReal, NULL, sizeof(qopt_Real), NULL };

    const char *str1 = NULL,
               *str2 = NULL;
    qopt_Int bool1 = 0;
    

    qoptElem qopt_tab_opt[] = {
        { QOPT_NEXT, NULL,      qOptArray,  0,  { .arr = &qopt_arr_cplx } },
        { QOPT_NEXT, NULL,      qOptArray,  0,  { .arr = &qopt_arr_str } },
        { QOPT_KEY, "struct1",  qOptStruct, 0,  { .tab = qopt_tab_struct1 } },
        { QOPT_KEY, "struct2",  qOptStruct, 0,  { .tab = qopt_tab_struct2 } },
        { QOPT_KEY, "int_Nx4",  qOptArray,  0,  { .arr = &qopt_arr_int1 } },
        { QOPT_KEY, "int_2x3x4",qOptArray,  0,  { .arr = &qopt_arr_int2 } },
        { QOPT_KEY, "real_Nx2", qOptArray,  0,  { .arr = &qopt_arr_real } },
        { QOPT_KEY, "str1",     qString,    0,  { .s = &str1 } },
        { QOPT_KEY, "str2",     qString,    QOPT_OPT,   { .s = &str2 } },
        { QOPT_KEY, "bool1",    qOptBool,   QOPT_OPT,   { .i = &bool1 } },
        { QOPT_END }
    } ;
    
    qoptElem stack_elem[] = {
        { QOPT_NEXT, NULL, qOptInt,     0,  { .i = &i1 } },
        { QOPT_NEXT, NULL, qOptInt,     0,  { .i = &i2 } },
        { QOPT_NEXT, NULL, qReal,       0,  { .r = &r3 } },
        { QOPT_NEXT, NULL, qComplex,    0,  { .c = &c4 } },
        { QOPT_NEXT, NULL, qOptArray,   0,  { .arr = &qopt_arr_recx } },
        { QOPT_NEXT, NULL, qOptStruct,  QOPT_OPT, { .tab = qopt_tab_opt } },
        { QOPT_END }
    };

    int status = 0;
    
    printf0("START lua_gettop(L) = %d\n", lua_gettop(L));
    status = qopt_read_stack(L, stack_elem);
    printf0("qopt read status = %d\n", status);
    printf0("STOP lua_gettop(L) = %d\n", lua_gettop(L));
    /* print parsed values */
    printf0("i1=%d  i2=%d  r3=%e  c4=(%e, %e)\n", 
            (int)i1, (int)i2, (double)r3, creal(c4), cimag(c4));
    printf0("struct1 = {%d, %d, %e, %d}\n", 
            s1i1, s1i2, s1r3, s1i4);
    printf0("struct2 = {%d, %e, %e, %e}\n",
            s2i1, s2r2, s2r3, s2r4);
    
    printf0("str1='%s'  str2='%s'  bool1=%d\n", str1, str2, (int)bool1);

    print_dim(2, real_dim, "real");
    print_arr_general(2, real_dim, qopt_arr_real.buf, sizeof(qopt_Real), print_qoptreal_, "real");

    print_dim(2, int1_dim, "int1");
    print_arr_general(2, int1_dim, qopt_arr_int1.buf, sizeof(qopt_Int), print_qoptint_, "int1");
    
    print_dim(3, int2_dim, "int2");
    print_arr_general(3, int2_dim, int2, sizeof(qopt_Int), print_qoptint_, "int2");

    print_dim(3, cplx_dim, "cplx");
    print_arr_general(3, cplx_dim, qopt_arr_cplx.buf, sizeof(qopt_Complex), print_qoptcomplex_, "cplx");
    
    print_dim(2, recx_dim, "recx");
    print_arr_general(2, recx_dim, qopt_arr_recx.buf, sizeof(recx), print_recx_, "recx");
    
    print_dim(2, str_dim, "str");
    print_arr_general(2, str_dim, str_arr_buf, sizeof(char *), print_str_, "str");

    /* freeing stuff */
    if (NULL != qopt_arr_real.buf) free(qopt_arr_real.buf);
    if (NULL != qopt_arr_int1.buf) free(qopt_arr_int1.buf);
    if (NULL != qopt_arr_cplx.buf) free(qopt_arr_cplx.buf);
    if (NULL != qopt_arr_recx.buf) free(qopt_arr_recx.buf);
  
    return 0;
}
