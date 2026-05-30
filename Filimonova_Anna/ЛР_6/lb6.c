#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void mantpor(float x, float *mant, double *por) {
    int st;
    if (x == 0.0f) {
        *mant = 0.0f;
        *por = 1.0;
        return;
    }
    st = (int)floor(log10(fabs(x)));
    *por = pow(10.0, st);
    *mant = x / (*por);
}

float round_mant(float x, int n) {
    if (x == 0.0f) return 0.0f;
    float xabs = fabs(x);
    int signx = (x >= 0) ? 1 : -1;
    float mant;
    double por;
    mantpor(xabs, &mant, &por);
    float dpown = pow(10.0f, n);
    float xr = signx * (floor(mant * dpown + 0.5f) / dpown) * por;
    return xr;
}

int miter(double** a, double* b, double* x, int n,
                  double eps, int max_iter, int* iter_count) {
    if (n < 2) return 2;

    double normA = 0.0;
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) sum += fabs(a[i][j]);
        if (sum > normA) normA = sum;
    }
    if (normA >= 1.0) return 1;

    for (int i = 0; i < n; i++) x[i] = b[i];

    double* x_old = (double*)malloc(n * sizeof(double));
    if (x_old == NULL) return 1;

    *iter_count = 0;
    double diff_norm, x_norm;

    do {
        for (int i = 0; i < n; i++) x_old[i] = x[i];

        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) sum += a[i][j] * x_old[j];
            x[i] = b[i] + sum;
        }

        (*iter_count)++;

        double diff2 = 0.0, norm2 = 0.0;
        for (int i = 0; i < n; i++) {
            double d = x[i] - x_old[i];
            diff2 += d * d;
            norm2 += x[i] * x[i];
        }
        diff_norm = sqrt(diff2);
        x_norm = sqrt(norm2);

        if (*iter_count > max_iter) {
            free(x_old);
            return 1;
        }
    } while (diff_norm / x_norm > eps);

    free(x_old);
    return 0;
}

void test_large_system() {
    int n = 10000;
    double eps = 1e-3;
    int max_iter = 100000;

    printf("\nТестирование большой системы: n = %d, eps = %g\n", n, eps);

    double** C = (double**)malloc(n * sizeof(double*));
    double** A = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        C[i] = (double*)calloc(n, sizeof(double));
        A[i] = (double*)calloc(n, sizeof(double));
    }
    double* d = (double*)malloc(n * sizeof(double));
    double* b = (double*)malloc(n * sizeof(double));
    double* x_true = (double*)malloc(n * sizeof(double));
    double* x = (double*)malloc(n * sizeof(double));

    srand(time(NULL));
    for (int i = 0; i < n; i++) {
        double sum_off = 0.0;
        for (int j = 0; j < n; j++) {
            if (i != j) {
                C[i][j] = (double)rand() / RAND_MAX * 2.0 - 1.0;
                sum_off += fabs(C[i][j]);
            }
        }
        C[i][i] = sum_off + 2.0 + (double)rand() / RAND_MAX;
    }

    for (int i = 0; i < n; i++)
        x_true[i] = (double)rand() / RAND_MAX * 20.0 - 10.0;

    for (int i = 0; i < n; i++) {
        d[i] = 0.0;
        for (int j = 0; j < n; j++)
            d[i] += C[i][j] * x_true[j];
        b[i] = d[i] / C[i][i];
        for (int j = 0; j < n; j++)
            A[i][j] = (i == j) ? 0.0 : -C[i][j] / C[i][i];
    }

    double normA = 0.0;
    for (int i = 0; i < n; i++) {
        double sum_row = 0.0;
        for (int j = 0; j < n; j++)
            sum_row += fabs(A[i][j]);
        if (sum_row > normA) normA = sum_row;
    }
    printf("||A||_inf = %.6f\n", normA);

    clock_t start = clock();
    int iter_count;
    int err = miter(A, b, x, n, eps, max_iter, &iter_count);
    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;

    if (err == 0) {
        double err_rel = 0.0, norm_true = 0.0;
        for (int i = 0; i < n; i++) {
            double diff = x[i] - x_true[i];
            err_rel += diff * diff;
            norm_true += x_true[i] * x_true[i];
        }
        err_rel = sqrt(err_rel) / sqrt(norm_true);
        printf("Сошелся за %d итераций\n", iter_count);
        printf("Время: %.3f сек\n", elapsed);
        printf("Относительная погрешность: %g\n", err_rel);
    } else {
        printf("Не сошелся (код ошибки %d)\n", err);
    }

    for (int i = 0; i < n; i++) {
        free(C[i]);
        free(A[i]);
    }
    free(C); free(A);
    free(d); free(b); free(x_true); free(x);
}

int main() {
    double A_small[3][3] = {
        { 0.0, -0.06,  0.02},
        {-0.03,  0.0,   0.05},
        {-0.01,  0.02,  0.0 }
    };
    double b_small[3] = {2.0, 3.0, 5.0};
    double x_small[3];
    int n_small = 3;
    int iter;

    double* A_rows[3] = {A_small[0], A_small[1], A_small[2]};

    FILE *fout = fopen("result.txt", "w");
    if (!fout) {
        perror("Cannot open output file");
        return 1;
    }

    printf("eps\t\titerations\n");
    fprintf(fout, "# eps  iterations\n");

    for (int k = 1; k <= 6; k++) {
        double eps = pow(10.0, -k);
        int err = miter(A_rows, b_small, x_small, n_small, eps, 10000, &iter);
        if (err != 0) {
            printf("%e\terror %d\n", eps, err);
            fprintf(fout, "%e\t-1\n", eps);
        } else {
            printf("%e\t%d\n", eps, iter);
            fprintf(fout, "%e\t%d\n", eps, iter);
        }
    }

    fclose(fout);

    int err_final = miter(A_rows, b_small, x_small, n_small, 1e-6, 10000, &iter);
    if (err_final == 0) {
        printf("\nSolution at eps=1e-6:\n");
        for (int i = 0; i < n_small; i++)
            printf("x[%d] = %.6f\n", i, x_small[i]);
    }

    srand(time(NULL));
    test_large_system();

    return 0;
}