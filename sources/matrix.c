#include <qlua.h>                                                    /* DEPS */
#include <matrix.h>                                                  /* DEPS */
#include <math.h>

/* symetric real matrix eigenvalues and eigenvectors */
#define a(M,i,j)   ((M)[((i)+(j)*n)])

static void
swap_mx(int i, int j, int k, int n, double *t)
{
    double x;
    int m;
    
    if (i == j)
        return;
    x = a(t,i,i); a(t,i,i) = a(t,j,j); a(t,j,j) = x;
    for (m = k; m < i; m++) {
        x = a(t,i,m); a(t,i,m) = a(t,j,m); a(t,j,m) = x;
    }
    for (m = i + 1; m < j; m++) {
        x = a(t,m,i); a(t,m,i) = a(t,j,m); a(t,j,m) = x;
    }
    for (m = j + 1; m < n; m++) {
        x = a(t,m,i); a(t,m,i) = a(t,m,j); a(t,m,j) = x;
    }
}

static void
swap_u(int i, int j, int n, double *u)
{
    double x;
    int m;

    for (m = 0; m < n; m++) {
        x = a(u,m,i); a(u,m,i) = a(u,m,j); a(u,m,j) = x;
    }
}

static void
tridiag(int n, double *d, double *e, double *t, double *u)
{
    int i;
    double K;

    /* order the lower triange of t in decreasing vector norm, adjust u */
    /*  compute column norms on t */
    for (i = 0; i < n; i++) {
        double v = fabs(a(t,i,i));
        int j;
        for (j = 0; j < i; j++) {
            v = hypot(v, a(t,i,j));
        }
        for (j = i + 1; j < n; j++) {
            v = hypot(v, a(t,j,i));
        }
        d[i] = v;
    }
    /*  reorder t and u */
    for (i = 0; i < n; i++) {
        double v = d[i];
        int j = i;
        int k;
        for (k = i + 1; k < n; k++) {
            if (d[k] > v) {
                j = k;
                v = d[k];
            }
        }
        d[j] = d[i];
        d[i] = v;
        swap_mx(i, j, 0, n, t);
        swap_u(i, j, n, u);
    }
    
    /* loop of householder transformations */
    for (i = 0; i < n - 1; i++) {
        double iv;
        int ii = i + 1;
        int im, jm;

        /* find max element in [ii..n-1,i] */
        iv = fabs(a(t,ii,i));
        for (im = ii, jm = ii + 1; jm < n; jm++) {
            double jv = fabs(a(t,jm,i));
            if (jv > iv) {
                im = jm;
                iv = jv;
            }
        }
        if (iv == 0) {
            e[i] = 0;
            d[i] = a(t,i,i);
            goto no_transform;
        }
        /* move max element to [ii,i] */
        swap_mx(ii,im, i, n, t);
        swap_u(ii,im, n, u);

        /* make t[ii,i] real */
        {
            double x = a(t,ii,i) / iv;
            int j;

            a(t,ii,i) = iv;
            for (j = ii + 1; j < n; j++) {
                double v = a(t,j,ii);
                a(t,j,ii) = v * x;
            }
            for (j = 0; j < n; j++) {
                double v = a(u,j,ii);
                a(u,j,ii) = v * x;
            }
        }
        /* find the householder transform and save d[i] and e[i] */
        d[i] = a(t,i,i);
        {
            double ar = a(t,ii,i);
            double dx = ar;
            int j;

            for (j = ii+1; j < n; j++) {
                dx = hypot(dx, a(t,j,i));
            }

            if (dx == ar) {
                e[i] = ar;
                goto no_transform;
            }
            e[i] = dx;
            a(t,ii,i) = (ar - dx);
            dx = 1 / sqrt(dx * (dx - ar));
            a(t,ii,i) *= dx;

            for (j = ii + 1; j < n; j++) {
                a(t,j,i) *= dx;
            }
        }
        /* multiply u by the transform */
        {
            int j;
            for (j = 0; j < n; j++) {
                double x = 0;
                int k;
                for (k = ii; k < n; k++) {
                    x += a(u,j,k) * a(t,k,i);
                }
                for (k = ii; k < n; k++) {
                    double v = a(t,k,i);
                    a(u,j,k) -= x * v;
                }
            }
        }
        /* apply the transform to t */
        /*  compute p in t[i,ii..] */
        {
            int j;
            for (j = ii; j < n; j++) {
                double x = 0;
                int k;

                for (k = ii; k < j; k++) {
                    double v = a(t,j,k);
                    x += v * a(t,k,i);
                }
                x += a(t,j,j) * a(t,j,i);
                for (k = j + 1; k < n; k++) {
                    double v = a(t,k,j);
                    x += v * a(t,k,i);
                }
                a(t,i,j) = x;
            }
        }
        /*  compute K */
        {
            int j;
            for (K = 0, j = ii; j < n; j++) {
                K += a(t,j,i) * a(t,i,j);
            }
            K = 0.5 * K;
        }
        /*  compute q in t[i,ii..] */
        {
            int j;
            for (j = ii; j < n; j++) {
                a(t,i,j) -= K * a(t,j,i);
            }
        }
        /*  compute A - uq* - qu* for L/ii */
        {
            int j;
            for (j = ii; j < n; j++) {
                double ujr = a(t,j,i);
                double qjr = a(t,i,j);
                int k;
                a(t,j,j) -= 2 * ujr * qjr;
                for (k = j + 1; k < n; k++) {
                    double ukr = a(t,k,i);
                    double qkr = a(t,i,k);
                    a(t,k,j) -= ujr*qkr + ukr*qjr;
                }
            }
        }
    no_transform:
        ;
    }
    /* handle the last row */
    d[n-1] = a(t, n-1, n-1);
}

static double
eigen2(double a, double b, double ab)
{
    double d = 0.5 * (a - b);

    if (d > 0) {
        return b - ab * (ab / (hypot(d, ab) + d));
    } else if (d < 0) {
        return b + ab * (ab / (hypot(d, ab) - d));
    } else {
        return b - fabs(ab);
    }
}

static void
givens(double *c, double *s, double a, double b)
{
    if (b == 0) {
        *c = 1;
        *s = 0;
    } else if (fabs(b) > fabs(a)) {
        double t = - a / b;
        double vs = 1/sqrt(1 + t * t);
        *c = t * vs;
        *s = vs;
    } else {
        double t = - b / a;
        double vc = 1/sqrt(1 + t * t);
        *c = vc;
        *s = t * vc;
    }
}

static void
rotate(int n, double *u, double c, double s, int j, int k)
{
    int i;

    for (i = 0; i < n; i++) {
        double ujr = a(u,i,j);
        double ukr = a(u,i,k);
        a(u,i,j) = c * ujr - s * ukr;
        a(u,i,k) = c * ukr + s * ujr;
    }
}

static void
qr_all(int n, double *d, double *e, double *u)
{
    int j, k = n - 1;

    for (;;) {
        double dx, mu, x, z, ap = 0, bp, aq = 0;
        double c, s;

        for (; k > 0; k--) {
            dx = fabs(d[k]) + fabs(d[k-1]);
            if (dx == dx + fabs(e[k-1]))
                e[k-1] = 0;
            else
                break;
        }
        if (k == 0)
            return;
        for (j = k - 1; j > 0; j--) {
            dx = fabs(d[j]) + fabs(d[j-1]);
            if (dx == dx + fabs(e[j-1])) {
                e[j-1] = 0;
                break;
            }
        }
        mu = eigen2(d[k - 1], d[k], e[k - 1]);
        ap = d[j];
        x = ap - mu;
        bp = z = e[j];
        aq = d[j + 1];
        if (k - j == 1) {
            double pcms, pspc, qsmc, qcps, ap1, bp1, aq1;
            givens(&c, &s, x, z);
            rotate(n, u, c, s, j, k);
            pcms = c * ap - s * bp;
            pspc = s * ap + c * bp;
            qsmc = s * aq - c * bp;
            qcps = c * aq + s * bp;
            ap1 = c * pcms + s * qsmc;
            bp1 = c * pspc - s * qcps;
            aq1 = s * pspc + c * qcps;
            d[j] = ap1;
            d[j + 1] = aq1;
            e[j] = bp1;
        } else {
            int m;
            double ak = 0, bk = 0, zk = 0, bq = e[j + 1];
            double bk1, pcms, qsmc, pspc, qcps;
            double ap1 = 0, bp1 = 0, zp1, aq1 = 0, bq1;

            for (m = j; m < k; m++) {
                givens(&c, &s, x, z);
                rotate(n, u, c, s, m, m + 1);
                bk1 = c * bk - s * zk;
                pcms = c * ap - s * bp;
                pspc = s * ap + c * bp;
                qsmc = s * aq - c * bp;
                qcps = c * aq + s * bp;
                ap1 = c * pcms + s * qsmc;
                bp1 = c * pspc - s * qcps;
                zp1 = - s * bq;
                aq1 = s * pspc + c * qcps;
                bq1 = c * bq;
                d[m] = ap1;
                if (m > j) e[m - 1] = bk1;
                if (m < k - 1) e[m + 1] = bq1;
                x = bp1;
                z = zp1;
                ak = ap1;
                bk = bp1;
                zk = zp1;
                ap = aq1;
                bp = bq1;
                if (m < k - 1) aq = d[m + 2];
                if (m < k - 2) bq = e[m + 2];
            }
            d[k] = aq1;
            e[k - 1] = bp1;
        }
    }
}

static void
order_lambdas(int n, double *d, double *u)
{
    int i;

    for (i = 0; i < n - 1; i++) {
        double x;
        int jv = i;
        int j;
        for (j = i + 1; j < n; j++) {
            if (d[j] < d[jv])
                jv = j;
        }
        x = d[jv]; d[jv] = d[i]; d[i] = x;
        swap_u(i, jv, n, u);
    }
}

void /* symmetic case */
matrix_reigenvec(lua_State *L,
                 int n, const double *data, double *lambda, double *u)
{
    double *t = qlua_malloc(L, n * n * sizeof (double));
    double *e = qlua_malloc(L, (n - 1) * sizeof (double));
    int i, j;

    for (i = 0; i < n; i++) {
        a(u, i, i) = 1.0;
        a(t, i, i) = a(data, i, i);
        for (j = 0; j < i; j++) {
            a(t, i, j) = a(data, i, j);
            a(t, j, i) = a(data, i, j);
            a(u, i, j) = 0;
            a(u, j, i) = 0;
        }
    }

    tridiag(n, lambda, e, t, u);
    qr_all(n, lambda, e, u);
    order_lambdas(n, lambda, u);
        
    qlua_free(L, t);
    qlua_free(L, e);
}

/* real matrix inversion */
int
matrix_rinverse(int n, double *M, double *N)
{
#define U(m,i,j) ((m)[(i)*(n)+(j)])
    int i, j, k;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            U(N,i,j) = 0;
        U(N,i,i) = 1;
    }
    
    for (i = 0; i < n; i++) {
        double m;
        for (k = i, m = fabs(U(M,i,i)), j = i + 1; j < n; j++) {
            double x = fabs(U(M,j,i));
            if (x > m) {
                m = x;
                k = j;
            }
        }
        if (m == 0)
            return 1;
        if (k != i) {
            for (j = i; j < n; j++) {
                double t = U(M,k,j);
                U(M,k,j) = U(M,i,j);
                U(M,i,j) = t;
            }
            for (j = 0; j < n; j++) {
                double t = U(N,k,j);
                U(N,k,j) = U(N,i,j);
                U(N,i,j) = t;
            }
        }
        m = 1.0 / U(M,i,i);
        for (j = i + 1; j < n; j++) {
            double a = m * U(M,j,i);
            for (k = i + 1; k < n; k++)
                U(M,j,k) = U(M,j,k) - a * U(M,i,k);
            for (k = 0; k < n; k++)
                U(N,j,k) = U(N,j,k) - a * U(N,i,k);
        }
    }

    for (i = n; i--;) {
        double m = 1.0 / U(M,i,i);
        for (k = 0; k < n; k++)
            U(N,i,k) *= m;
        for (j = 0; j < i; j++) {
            double a = U(M,j,i);
            for (k = 0; k < n; k++) {
                U(N,j,k) = U(N,j,k) - a * U(N,i,k);
            }
        }
    }

    return 0;
#undef U
}
