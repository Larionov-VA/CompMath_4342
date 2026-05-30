#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define nmax 10000
#define EPS 1e-12

int gauss(double **a, double *b, double *x, int n, double eps) {
    if (n < 2) return 2;
    if (n > nmax) return 3;
    int i, j, k, imax;
    double amax, t, mi, sum;
    for (k = 0; k < n-1; k++) {
        imax = k;
        amax = fabs(a[k][k]);
        for (i = k+1; i < n; i++)
            if (fabs(a[i][k]) > amax) {
                amax = fabs(a[i][k]);
                imax = i;
            }
        if (amax < eps) return 1;
        if (imax != k) {
            for (j = k; j < n; j++) {
                t = a[imax][j];
                a[imax][j] = a[k][j];
                a[k][j] = t;
            }
            t = b[imax];
            b[imax] = b[k];
            b[k] = t;
        }
        for (i = k+1; i < n; i++) {
            mi = a[i][k] / a[k][k];
            a[i][k] = 0.0;
            for (j = k+1; j < n; j++)
                a[i][j] -= mi * a[k][j];
            b[i] -= mi * b[k];
        }
    }
    if (fabs(a[n-1][n-1]) < eps) return 1;
    for (i = n-1; i >= 0; i--) {
        sum = 0.0;
        for (j = i+1; j < n; j++) sum += a[i][j] * x[j];
        x[i] = (b[i] - sum) / a[i][i];
    }
    return 0;
}

double residual(double **A, double *b, double *x, int n) {
    double res = 0.0;
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) sum += A[i][j] * x[j];
        double diff = sum - b[i];
        res += diff * diff;
    }
    return sqrt(res);
}

double **load_matrix_txt(const char *filename, int *n) {
    FILE *f = fopen(filename, "r");
    if (!f) return NULL;
    fscanf(f, "%d", n);
    double **A = malloc(*n * sizeof(double*));
    for (int i = 0; i < *n; i++) {
        A[i] = malloc(*n * sizeof(double));
        for (int j = 0; j < *n; j++)
            fscanf(f, "%lf", &A[i][j]);
    }
    fclose(f);
    return A;
}

void benchmark_random(int n) {
    char fname[64];
    sprintf(fname, "random_%d.txt", n);
    int n_loaded;
    double **A = load_matrix_txt(fname, &n_loaded);
    if (!A || n_loaded != n) return;

    double **A_copy = malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        A_copy[i] = malloc(n * sizeof(double));
        for (int j = 0; j < n; j++) A_copy[i][j] = A[i][j];
    }

    double *x_true = malloc(n * sizeof(double));
    double *b = malloc(n * sizeof(double));
    double *b_copy = malloc(n * sizeof(double));
    for (int i = 0; i < n; i++)
        x_true[i] = (double)rand() / RAND_MAX * 20.0 - 10.0;
    for (int i = 0; i < n; i++) {
        b[i] = 0.0;
        for (int j = 0; j < n; j++) b[i] += A[i][j] * x_true[j];
        b_copy[i] = b[i];
    }

    double *x = malloc(n * sizeof(double));
    clock_t start = clock();
    int err = gauss(A, b, x, n, EPS);
    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;

    if (err == 0) {
        double rel_err = 0.0, norm_true = 0.0;
        for (int i = 0; i < n; i++) {
            double d = x[i] - x_true[i];
            rel_err += d*d;
            norm_true += x_true[i]*x_true[i];
        }
        rel_err = sqrt(rel_err) / sqrt(norm_true);
        double res = residual(A_copy, b_copy, x, n);
        printf("R %d %.3f %.6e %.6e\n", n, elapsed, rel_err, res);
    } else {
        printf("R %d -1.0 -1.0 -1.0\n", n);
    }

    for (int i = 0; i < n; i++) { free(A[i]); free(A_copy[i]); }
    free(A); free(A_copy); free(b); free(b_copy); free(x); free(x_true);
}

void benchmark_hilbert(int n) {
    char fname[64];
    sprintf(fname, "hilbert_%d.txt", n);
    int n_loaded;
    double **H = load_matrix_txt(fname, &n_loaded);
    if (!H || n_loaded != n) return;

    double **H_copy = malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        H_copy[i] = malloc(n * sizeof(double));
        for (int j = 0; j < n; j++) H_copy[i][j] = H[i][j];
    }

    double *x_true = malloc(n * sizeof(double));
    double *b = malloc(n * sizeof(double));
    double *b_copy = malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) x_true[i] = 1.0;
    for (int i = 0; i < n; i++) {
        b[i] = 0.0;
        for (int j = 0; j < n; j++) b[i] += H[i][j] * x_true[j];
        b_copy[i] = b[i];
    }

    double *x = malloc(n * sizeof(double));
    int err = gauss(H, b, x, n, EPS);
    if (err == 0) {
        double rel_err = 0.0, norm_true = 0.0;
        for (int i = 0; i < n; i++) {
            double d = x[i] - x_true[i];
            rel_err += d*d;
            norm_true += x_true[i]*x_true[i];
        }
        rel_err = sqrt(rel_err) / sqrt(norm_true);
        double res = residual(H_copy, b_copy, x, n);
        printf("H %d %.6e %.6e\n", n, rel_err, res);
    } else {
        printf("H %d -1.0 -1.0\n", n);
    }

    for (int i = 0; i < n; i++) { free(H[i]); free(H_copy[i]); }
    free(H); free(H_copy); free(b); free(b_copy); free(x); free(x_true);
}

void test_basic() {
    double **A = malloc(2 * sizeof(double*));
    double **A_copy = malloc(2 * sizeof(double*));
    for (int i = 0; i < 2; i++) {
        A[i] = malloc(2 * sizeof(double));
        A_copy[i] = malloc(2 * sizeof(double));
    }
    double b[2], b_copy[2], x[2];
    A[0][0] = 2; A[0][1] = 1; b[0] = 3;
    A[1][0] = 1; A[1][1] = 1; b[1] = 2;
    for (int i = 0; i < 2; i++) {
        b_copy[i] = b[i];
        for (int j = 0; j < 2; j++) A_copy[i][j] = A[i][j];
    }
    int err = gauss(A, b, x, 2, EPS);
    if (err == 0) {
        double res = residual(A_copy, b_copy, x, 2);
        printf("2x2: (%.12f, %.12f) residual %.6e\n", x[0], x[1], res);
    } else {
        printf("2x2: error %d\n", err);
    }
    for (int i = 0; i < 2; i++) { free(A[i]); free(A_copy[i]); }
    free(A); free(A_copy);
}

void test_singular() {
    double **A = malloc(2 * sizeof(double*));
    for (int i = 0; i < 2; i++) A[i] = malloc(2 * sizeof(double));
    double b[2], x[2];
    A[0][0] = 1; A[0][1] = 2; b[0] = 3;
    A[1][0] = 2; A[1][1] = 4; b[1] = 6;
    int err = gauss(A, b, x, 2, EPS);
    printf("SINGULAR: %s\n", err == 1 ? "detected" : (err == 0 ? "not detected" : "error"));
    for (int i = 0; i < 2; i++) free(A[i]);
    free(A);
}

int main() {
    srand(time(NULL));
    test_basic();
    test_singular();

    int sizes[] = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
    for (int i = 0; i < 10; i++) benchmark_random(sizes[i]);

    int hsizes[] = {2,3,4,5,6,7,8,9,10,11};
    for (int i = 0; i < 9; i++) benchmark_hilbert(hsizes[i]);

    return 0;
}