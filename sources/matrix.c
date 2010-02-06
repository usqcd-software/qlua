#include <qlua.h>                                                    /* DEPS */
#include <matrix.h>                                                  /* DEPS */
#include <math.h>

/* symetric real matrix eigenvalues and eigenvectors */
#define a(M,i,j)   ((M)[((i)+(j)*n)])

static void
r_swap_mx(int i, int j, int k, int n, double *t)
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
r_swap_u(int i, int j, int n, double *u)
{
    double x;
    int m;

    for (m = 0; m < n; m++) {
        x = a(u,m,i); a(u,m,i) = a(u,m,j); a(u,m,j) = x;
    }
}

static void
r_tridiag(int n, double *d, double *e, double *t, double *u)
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
        r_swap_mx(i, j, 0, n, t);
        r_swap_u(i, j, n, u);
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
        r_swap_mx(ii,im, i, n, t);
        r_swap_u(ii,im, n, u);

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
r_rotate(int n, double *u, double c, double s, int j, int k)
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
r_qr_all(int n, double *d, double *e, double *u)
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
            r_rotate(n, u, c, s, j, k);
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
                r_rotate(n, u, c, s, m, m + 1);
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
r_order_lambdas(int n, double *d, double *u)
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
        r_swap_u(i, jv, n, u);
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

    r_tridiag(n, lambda, e, t, u);
    r_qr_all(n, lambda, e, u);
    r_order_lambdas(n, lambda, u);
        
    qlua_free(L, t);
    qlua_free(L, e);
}
#undef a

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

/* hermitian complex matrix eigenvalues and eigenvectors */
#define a(M,i,j)   (&((M)[2*((i)+(j)*n)]))

static void
c_swap_mx(int i, int j, int k, int n, double *t)
{
    double x;
    int m;
    
    if (i == j)
        return;
    a(t,j,i)[1] = - a(t,j,i)[1];
    x = a(t,i,i)[0]; a(t,i,i)[0] = a(t,j,j)[0]; a(t,j,j)[0] = x;
    for (m = k; m < i; m++) {
        x = a(t,i,m)[0]; a(t,i,m)[0] = a(t,j,m)[0]; a(t,j,m)[0] = x;
        x = a(t,i,m)[1]; a(t,i,m)[1] = a(t,j,m)[1]; a(t,j,m)[1] = x;
    }
    for (m = i + 1; m < j; m++) {
        x = a(t,m,i)[0]; a(t,m,i)[0] = a(t,j,m)[0]; a(t,j,m)[0] = x;
        x = a(t,m,i)[1]; a(t,m,i)[1] = -a(t,j,m)[1]; a(t,j,m)[1] = -x;
    }
    for (m = j + 1; m < n; m++) {
        x = a(t,m,i)[0]; a(t,m,i)[0] = a(t,m,j)[0]; a(t,m,j)[0] = x;
        x = a(t,m,i)[1]; a(t,m,i)[1] = a(t,m,j)[1]; a(t,m,j)[1] = x;
    }
}

static void
c_swap_u(int i, int j, int n, double *u)
{
    double x;
    int m;

    for (m = 0; m < n; m++) {
        x = a(u,m,i)[0]; a(u,m,i)[0] = a(u,m,j)[0]; a(u,m,j)[0] = x;
        x = a(u,m,i)[1]; a(u,m,i)[1] = a(u,m,j)[1]; a(u,m,j)[1] = x;
    }
}

static void
c_tridiag(int n, double *d, double *e, double *t, double *u)
{
    int i;
    double K;

    /* order the lower triange of t in decreasing vector norm, adjust u */
    /*  compute column norms on t */
    for (i = 0; i < n; i++) {
        double v = fabs(a(t,i,i)[0]);
        int j;
        for (j = 0; j < i; j++) {
            v = hypot(v, a(t,i,j)[0]);
            v = hypot(v, a(t,i,j)[1]);
        }
        for (j = i + 1; j < n; j++) {
            v = hypot(v, a(t,j,i)[0]);
            v = hypot(v, a(t,j,i)[1]);
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
        c_swap_mx(i, j, 0, n, t);
        c_swap_u(i, j, n, u);
    }
    
    /* loop of householder transformations */
    for (i = 0; i < n - 1; i++) {
        double iv;
        int ii = i + 1;
        int im, jm;

        /* find max element in [ii..n-1,i] */
        iv = hypot(a(t,ii,i)[0], a(t,ii,i)[1]);
        for (im = ii, jm = ii + 1; jm < n; jm++) {
            double jv = hypot(a(t,jm,i)[0], a(t,jm,i)[1]);
            if (jv > iv) {
                im = jm;
                iv = jv;
            }
        }
        if (iv == 0) {
            e[i] = 0;
            d[i] = a(t,i,i)[0];
            goto no_transform;
        }
        /* move max element to [ii,i] */
        c_swap_mx(ii,im, i, n, t);
        c_swap_u(ii,im, n, u);

        /* make t[ii,i] real */
        {
            double x = a(t,ii,i)[0] / iv;
            double y = a(t,ii,i)[1] / iv;
            int j;

            a(t,ii,i)[0] = iv;
            a(t,ii,i)[1] = 0;
            for (j = ii + 1; j < n; j++) {
                double v = a(t,j,ii)[0];
                double w = a(t,j,ii)[1];
                a(t,j,ii)[0] = v * x - w * y;
                a(t,j,ii)[1] = v * y + w * x;
            }
            for (j = 0; j < n; j++) {
                double v = a(u,j,ii)[0];
                double w = a(u,j,ii)[1];
                a(u,j,ii)[0] = v * x - w * y;
                a(u,j,ii)[1] = v * y + w * x;
            }
        }
        /* find the householder transform and save d[i] and e[i] */
        d[i] = a(t,i,i)[0];
        {
            double ar = a(t,ii,i)[0];
            double dx = ar;
            int j;

            for (j = ii+1; j < n; j++) {
                dx = hypot(dx, a(t,j,i)[0]);
                dx = hypot(dx, a(t,j,i)[1]);
            }

            if (dx == ar) {
                e[i] = ar;
                goto no_transform;
            }
            e[i] = dx;
            a(t,ii,i)[0] = (ar - dx);
            dx = 1 / sqrt(dx * (dx - ar));
            a(t,ii,i)[0] *= dx;

            for (j = ii + 1; j < n; j++) {
                a(t,j,i)[0] *= dx;
                a(t,j,i)[1] *= dx;
            }
        }
        /* multiply u by the transform */
        {
            int j;
            for (j = 0; j < n; j++) {
                double x = 0;
                double y = 0;
                int k;
                for (k = ii; k < n; k++) {
                    x += a(u,j,k)[0] * a(t,k,i)[0] - a(u,j,k)[1] * a(t,k,i)[1];
                    y += a(u,j,k)[1] * a(t,k,i)[0] + a(u,j,k)[0] * a(t,k,i)[1];
                }
                for (k = ii; k < n; k++) {
                    double v = a(t,k,i)[0];
                    double w = a(t,k,i)[1];
                    a(u,j,k)[0] -= x * v + y * w;
                    a(u,j,k)[1] -= y * v - x * w;
                }
            }
        }
        /* apply the transform to t */
        /*  compute p in t[i,ii..] */
        {
            int j;
            for (j = ii; j < n; j++) {
                double x = 0;
                double y = 0;
                int k;

                for (k = ii; k < j; k++) {
                    double v = a(t,j,k)[0];
                    double w = a(t,j,k)[1];
                    x += v * a(t,k,i)[0] - w * a(t,k,i)[1];
                    y += v * a(t,k,i)[1] + w * a(t,k,i)[0];
                }
                x += a(t,j,j)[0] * a(t,j,i)[0];
                y += a(t,j,j)[0] * a(t,j,i)[1];
                for (k = j + 1; k < n; k++) {
                    double v = a(t,k,j)[0];
                    double w = a(t,k,j)[1];
                    x += v * a(t,k,i)[0] + w * a(t,k,i)[1];
                    y += v * a(t,k,i)[1] - w * a(t,k,i)[0];
                }
                a(t,i,j)[0] = x;
                a(t,i,j)[1] = y;
            }
        }
        /*  compute K */
        {
            int j;
            for (K = 0, j = ii; j < n; j++) {
                K += a(t,j,i)[0] * a(t,i,j)[0] + a(t,j,i)[1] * a(t,i,j)[1];
            }
            K = 0.5 * K;
        }
        /*  compute q in t[i,ii..] */
        {
            int j;
            for (j = ii; j < n; j++) {
                a(t,i,j)[0] -= K * a(t,j,i)[0];
                a(t,i,j)[1] -= K * a(t,j,i)[1];
            }
        }
        /*  compute A - uq* - qu* for L/ii */
        {
            int j;
            for (j = ii; j < n; j++) {
                double ujr = a(t,j,i)[0];
                double uji = a(t,j,i)[1];
                double qjr = a(t,i,j)[0];
                double qji = a(t,i,j)[1];
                int k;
                a(t,j,j)[0] -= 2 * (ujr * qjr + uji * qji);
                for (k = j + 1; k < n; k++) {
                    double ukr = a(t,k,i)[0];
                    double uki = a(t,k,i)[1];
                    double qkr = a(t,i,k)[0];
                    double qki = a(t,i,k)[1];
                    a(t,k,j)[0] -= ujr*qkr + uji*qki + ukr*qjr + uki*qji;
                    a(t,k,j)[1] -= ujr*qki - uji*qkr - ukr*qji + uki*qjr;
                }
            }
        }
    no_transform:
        ;
    }
    /* handle the last row */
    d[n-1] = a(t, n-1, n-1)[0];
}

static void
c_rotate(int n, double *u, double c, double s, int j, int k)
{
    int i;

    for (i = 0; i < n; i++) {
        double ujr = a(u,i,j)[0];
        double uji = a(u,i,j)[1];
        double ukr = a(u,i,k)[0];
        double uki = a(u,i,k)[1];
        a(u,i,j)[0] = c * ujr - s * ukr;
        a(u,i,j)[1] = c * uji - s * uki;
        a(u,i,k)[0] = c * ukr + s * ujr;
        a(u,i,k)[1] = c * uki + s * uji;
    }
}

static void
c_qr_all(int n, double *d, double *e, double *u)
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
            c_rotate(n, u, c, s, j, k);
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
                c_rotate(n, u, c, s, m, m + 1);
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
c_order_lambdas(int n, double *d, double *u)
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
        c_swap_u(i, jv, n, u);
    }
}

void /* hermitian case */
matrix_ceigenvec(lua_State *L,
                 int n, const double *data, double *lambda, double *u)
{
    double *t = qlua_malloc(L, 2 * n * n * sizeof (double));
    double *e = qlua_malloc(L, (n - 1) * sizeof (double));
    int i, j;

    for (i = 0; i < n; i++) {
        a(u, i, i)[0] = 1.0;
        a(u, i, i)[1] = 0.0;
        a(t, i, i)[0] = a(data, i, i)[0];
        a(t, i, i)[1] = 0;
        for (j = 0; j < i; j++) {
            a(t, i, j)[0] = a(data, i, j)[0];
            a(t, i, j)[1] = a(data, i, j)[1];
            a(t, j, i)[0] = a(data, i, j)[0];
            a(t, j, i)[1] = -a(data, i, j)[1];
            a(u, i, j)[0] = 0;
            a(u, i, j)[1] = 0;
            a(u, j, i)[0] = 0;
            a(u, j, i)[1] = 0;
        }
    }

    c_tridiag(n, lambda, e, t, u);
    c_qr_all(n, lambda, e, u);
    c_order_lambdas(n, lambda, u);
    
    qlua_free(L, t);
    qlua_free(L, e);
}
#undef a

/* complex matrix inversion */
int
matrix_cinverse(int n, double *M, double *N)
{
#define U(m,i,j) (&((m)[2*((i)*(n)+(j))]))
    int i, j, k;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            U(N,i,j)[0] = 0;
            U(N,i,j)[1] = 0;
        }
        U(N,i,i)[0] = 1;
        U(N,i,i)[1] = 0;
    }
    
    for (i = 0; i < n; i++) {
        double m, mr, mi;
        for (k = i, m = hypot(U(M,i,i)[0], U(M,i,i)[1]), j = i + 1;
             j < n; j++) {
            double x = hypot(U(M,j,i)[0], U(M,j,i)[1]);
            if (x > m) {
                m = x;
                k = j;
            }
        }
        if (m == 0)
            return 1;
        if (k != i) {
            for (j = i; j < n; j++) {
                double tr = U(M,k,j)[0];
                double ti = U(M,k,j)[1];
                U(M,k,j)[0] = U(M,i,j)[0];
                U(M,k,j)[1] = U(M,i,j)[1];
                U(M,i,j)[0] = tr;
                U(M,i,j)[1] = ti;
            }
            for (j = 0; j < n; j++) {
                double tr = U(N,k,j)[0];
                double ti = U(N,k,j)[1];
                U(N,k,j)[0] = U(N,i,j)[0];
                U(N,k,j)[1] = U(N,i,j)[1];
                U(N,i,j)[0] = tr;
                U(N,i,j)[1] = ti;
            }
        }
        m = 1.0 / hypot(U(M,i,i)[0], U(M,i,i)[1]);
        mr = U(M,i,i)[0] * m * m;
        mi = -U(M,i,i)[1] * m * m;
        for (j = i + 1; j < n; j++) {
            double ar = mr * U(M,j,i)[0] - mi * U(M,j,i)[1];
            double ai = mi * U(M,j,i)[0] + mr * U(M,j,i)[1];
            for (k = i + 1; k < n; k++) {
                U(M,j,k)[0] = U(M,j,k)[0] - ar * U(M,i,k)[0] + ai * U(M,i,k)[1];
                U(M,j,k)[1] = U(M,j,k)[1] - ai * U(M,i,k)[0] - ar * U(M,i,j)[1];
            }
            for (k = 0; k < n; k++) {
                U(N,j,k)[0] = U(N,j,k)[0] - ar * U(N,i,k)[0] + ai * U(N,i,k)[1];
                U(N,j,k)[1] = U(N,j,k)[1] - ai * U(N,i,k)[0] - ar * U(N,i,k)[1];
            }
        }
    }

    for (i = n; i--;) {
        double m = 1.0 / hypot(U(M,i,i)[0], U(M,i,i)[1]);
        double mr = U(M,i,i)[0] * m * m;
        double mi = -U(M,i,i)[1] * m * m;
        for (k = 0; k < n; k++) {
            double tr = U(N,i,k)[0] * mr - U(N,i,k)[1] * mi;
            double ti = U(N,i,k)[0] * mi + U(N,i,k)[1] * mr;
            U(N,i,k)[0] = tr;
            U(N,i,k)[1] = ti;
        }
        for (j = 0; j < i; j++) {
            double ar = U(M,j,i)[0];
            double ai = U(M,j,i)[1];
            for (k = 0; k < n; k++) {
                U(N,j,k)[0] = U(N,j,k)[0] - ar * U(N,i,k)[0] + ai * U(N,i,k)[1];
                U(N,j,k)[1] = U(N,j,k)[1] - ai * U(N,i,k)[0] - ar * U(N,i,k)[1];
            }
        }
    }

    return 0;
#undef U
}
