#ifndef MARK_E9EC2494_23C5_4621_898A_4414F4122DB2
#define MARK_E9EC2494_23C5_4621_898A_4414F4122DB2

/* real matrices */
int matrix_rinverse(int n, double *in, double *out);
void matrix_rqr(int n, int m, double *m_r, double *q);
/* symmetic case */
void matrix_reigenvec(lua_State *L, int n, double *u, double *lambda);

/* complex matrices */
int matrix_cinverse(int n, double *in, double *out);
void matrix_cqr(int n, int m, double *m_r, double *q);
/* hermitian case */
void matrix_ceigenvec(lua_State *L, int n, double *u, double *lambda);

#endif /* !defined(MARK_E9EC2494_23C5_4621_898A_4414F4122DB2) */
