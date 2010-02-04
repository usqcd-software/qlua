#ifndef MARK_E9EC2494_23C5_4621_898A_4414F4122DB2
#define MARK_E9EC2494_23C5_4621_898A_4414F4122DB2

void /* symmetic case */
matrix_reigenvec(lua_State *L,
                 int n, const double *data, double *lambda, double *u);
int matrix_rinverse(int n, double *in, double *out);

#endif /* !defined(MARK_E9EC2494_23C5_4621_898A_4414F4122DB2) */
