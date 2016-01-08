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
    test_parse(
      l_nev, l_ncv, l_maxiter, l_tol,
      { ["eigcg"] = {eigcg_vmax, eigcg_nev, eigcg_tol, eigcg_umax},
        ["cheb_accel"] = { l_poly_n, l_poly_a, l_poly_b, l_poly_x0 },
        ["coeff_list"]  = coeffs_list_Nx2,
        ["mom_list"]    = mom_list_Mx4,
        ["xy"] = xy_2x3x4,
        ["which"] = "LR", -- need the largest ev of polynomial T_n(A)
        ["arpack_logfile"] = arpack_logfile(cfg_key),
        ["inplace"] = true
        }
        abc_list)
    */
    qopt_Int l_nev = -1,
              l_ncv = -1,
              l_maxiter = -1;
    qopt_Real l_tol = -1.;


    /* variables to fill, with default values ; 
       while default values are not necessary for mandatory options, 
       the table itself will be mandatory, so better initialize */

    qopt_Int eigcg_vmax  = -1, 
              eigcg_nev   = -1, 
              eigcg_umax  = -1;
    qopt_Real eigcg_tol = -1.;
    qoptElem qopt_tab_eigcg[] = {
        { 1, NULL, qOptInt,   0, { .i = &eigcg_vmax } },
        { 2, NULL, qOptInt,   0, { .i = &eigcg_nev  } },
        { 3, NULL, qReal,  0, { .r = &eigcg_tol  } },
        { 4, NULL, qOptInt,   0, { .i = &eigcg_umax } },
        { QOPT_END }
    };
    
    qopt_Int l_n = -1;
    qopt_Real l_a  = DBL_MAX,
           l_b  = DBL_MAX,
           l_x0 = DBL_MAX;
    qoptElem qopt_tab_cheb_accel[] = {
        { 1, NULL, qOptInt,   0,   { .i = &l_n } },
        { 2, NULL, qReal,  0,   { .r = &l_a } },
        { 3, NULL, qReal,  0,   { .r = &l_b } },
        { 4, NULL, qReal,  QOPT_OPT,   { .r = &l_x0 } },
        { QOPT_END }
    };
    
    int coeffs_dim[]    = { -1, 2 },
        coeffs_dim_max[] = { 10, -1 }; 
    qoptNdarray qopt_arr_coeffs = { 2, coeffs_dim, coeffs_dim_max, 
                qReal, NULL, sizeof(qopt_Real), NULL };

    int mom_dim[] = {-1, 4},
        mom_dim_max[] = { 10, -1};
    qoptNdarray qopt_arr_mom = qopt_array_scalar(2, mom_dim, mom_dim_max, 
                qOptInt, NULL, NULL);

    int xy_dim[] = { 2, 3, 4 };
    qopt_Int xy[2][3][4];
    qoptNdarray qopt_arr_xy = qopt_array_scalar(3, xy_dim, NULL, 
                qOptInt, xy, NULL);

    const char *which = NULL;
    const char *arpack_logfile = NULL;
    qopt_Int inplace = 0;
    
    int abc_dim[]    = { -1, -1, -1 },
        abc_dim_max[] = { 10, 10, 10 };
    qoptNdarray qopt_arr_abc = qopt_array_scalar(3, abc_dim, abc_dim_max, 
                qComplex, NULL, NULL);

    int recx_dim[] = { -1, -1 };
    int recx_dim_max[] = { -1, -1 };
    qoptNdarray qopt_arr_recx = qopt_array_struct(2, recx_dim, recx_dim_max, 
                qopt_tab_recx, sizeof(recx), NULL, NULL);

    int str_dim[]     = { -1, -1 },
        str_dim_max[] = {  4,  4 };
    const char *str_arr_buf[16];
    qoptNdarray qopt_arr_str = qopt_array_scalar(2, str_dim, str_dim_max, 
                qString, str_arr_buf, NULL);


    qoptElem qopt_tab_opt[] = {
        { QOPT_KEY, "eigcg",       qOptStruct,  0, { .tab = qopt_tab_eigcg } },
        { QOPT_KEY, "cheb_accel",  qOptStruct,  0, { .tab = qopt_tab_cheb_accel } },
        { QOPT_KEY, "xy",          qOptArray,  0, { .arr = &qopt_arr_xy } },
        { QOPT_KEY, "mom_list",    qOptArray,  0, { .arr = &qopt_arr_mom } },
        { QOPT_KEY, "coeff_list",  qOptArray,  0, { .arr = &qopt_arr_coeffs } },
        { QOPT_KEY, "abc",         qOptArray,  0, { .arr = &qopt_arr_abc } },
        { QOPT_KEY, "which",       qString, 0, { .s = &which } },
        { QOPT_KEY, "arpack_logfile", qString, QOPT_OPT,
                    { .s = &arpack_logfile } },
        { QOPT_KEY, "inplace",     qOptBool, QOPT_OPT,
                    { .i = &inplace } },
        { QOPT_END }
    } ;
    
    qoptElem stack_elem[] = {
        { QOPT_NEXT, NULL, qOptArray, 0, { .arr = &qopt_arr_recx } },
        { QOPT_NEXT, NULL, qOptInt,   0, { .i = &l_nev } },
        { QOPT_NEXT, NULL, qOptInt,   0, { .i = &l_ncv } },
        { QOPT_NEXT, NULL, qOptInt,   0, { .i = &l_maxiter } },
        { QOPT_NEXT, NULL, qReal,  0, { .r = &l_tol } },
        { QOPT_NEXT, NULL, qOptArray, 0, { .arr = &qopt_arr_str } },
        { QOPT_NEXT, NULL, qOptStruct, QOPT_OPT, { .tab = qopt_tab_opt } },
        { QOPT_END }
    };

    int status = 0;
    
    printf0("START lua_gettop(L) = %d\n", lua_gettop(L));
    status = qopt_read_stack(L, stack_elem);
    printf0("qopt read status = %d\n", status);
    printf0("STOP lua_gettop(L) = %d\n", lua_gettop(L));
    /* print parsed values */
    printf0("l_nev=%d  l_ncv=%d  l_maxiter=%d  l_tol=%e\n", 
            l_nev, l_ncv, l_maxiter, l_tol);
    printf0("eigcg = {%d, %d, %e, %d}\n", 
            eigcg_vmax, eigcg_nev, eigcg_tol, eigcg_umax);
    printf0("cheb_accel = {%d, %e, %e, %e}\n",
            l_n, l_a, l_b, l_x0);
    
    printf0("which='%s'  arpack_logfile='%s'  inplace=%d\n", which, arpack_logfile, (int)inplace);

    print_dim(2, coeffs_dim, "coeffs");
    print_arr_general(2, coeffs_dim, qopt_arr_coeffs.buf, sizeof(qopt_Real), print_qoptreal_, "coeffs");

    print_dim(2, mom_dim, "mom");
    print_arr_general(2, mom_dim, qopt_arr_mom.buf, sizeof(qopt_Int), print_qoptint_, "mom");
    
    print_dim(3, xy_dim, "xy");
    print_arr_general(3, xy_dim, xy, sizeof(qopt_Int), print_qoptint_, "xy");

    print_dim(3, abc_dim, "abc");
    print_arr_general(3, abc_dim, qopt_arr_abc.buf, sizeof(qopt_Complex), print_qoptcomplex_, "abc");
    
    print_dim(2, recx_dim, "recx");
    print_arr_general(2, recx_dim, qopt_arr_recx.buf, sizeof(recx), print_recx_, "recx");
    
    print_dim(2, str_dim, "str");
    print_arr_general(2, str_dim, str_arr_buf, sizeof(char *), print_str_, "str");
  
    return 0;
}
