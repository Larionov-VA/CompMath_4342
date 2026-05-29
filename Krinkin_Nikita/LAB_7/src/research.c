#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#define ABS(num) (num < 0 ? -num : num)
#define M_TYPE float

typedef struct {
    M_TYPE **data;
    int n;
} matrix_t;

void print_table(matrix_t *matt) {
    if (matt) {
        for (int i=0; i<matt->n; i++) {
            for (int j=0; j<matt->n; j++) {
                if (sizeof(M_TYPE) == sizeof(float)) printf("%8f ", matt->data[i][j]);
                if (sizeof(M_TYPE) == sizeof(double)) printf("%8lf ", matt->data[i][j]);
            }
            // print b[i]
            if (sizeof(M_TYPE) == sizeof(float)) printf("| %8f ", matt->data[i][matt->n]);
            if (sizeof(M_TYPE) == sizeof(double)) printf("| %8lf ", matt->data[i][matt->n]);
            printf("\n");
        }
    }
}

matrix_t *create_empty(int n) {
    matrix_t *mat = malloc(sizeof(matrix_t));
    if (mat) {
        mat->n = n;
        mat->data = calloc(n, sizeof(M_TYPE **));
        for (int i=0; i<n; i++) {
            ((M_TYPE **) mat->data)[i] = (M_TYPE *) calloc(n + 1, sizeof(M_TYPE));
        }
    }
    return mat;
}

matrix_t *generate_random_sols(int n, M_TYPE *sols) {
    matrix_t *mat = create_empty(n);
    if (n) {
        // fill mat[i][j]
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                if (sizeof(M_TYPE) == sizeof(float)) {
                    mat->data[i][j] = (float) rand() / (float) rand();
                }
                if (sizeof(M_TYPE) == sizeof(double)) {
                    mat->data[i][j] = (float) rand() / (float) rand();
                }
            }
        }
        // generate sols
        for (int i=0; i<n; i++) {
            if (sizeof(M_TYPE) == sizeof(float)) {
                sols[i] = (float) rand() / (float) rand();
            }
            if (sizeof(M_TYPE) == sizeof(double)) {
                sols[i] = (float) rand() / (float) rand();
            }
        }
        // calculate b[i]
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                mat->data[i][n] += sols[j] * mat->data[i][j];

            }
        }
    }
    return mat;
}

matrix_t *generate_random(int n) {
    matrix_t *mat = create_empty(n);
    if (n) {
        for (int i=0; i<n; i++) {
            for (int j=0; j<n + 1; j++) {
                if (sizeof(M_TYPE) == sizeof(float)) {
                    mat->data[i][j] = (float) rand() / (float) rand();
                }
                if (sizeof(M_TYPE) == sizeof(double)) {
                    mat->data[i][j] = (float) rand() / (float) rand();
                }
            }
        }
    }
    return mat;
}

matrix_t *from_file(const char *path) {
    matrix_t *mat = 0;
    FILE *file = fopen(path, "r");
    if (file) {
        int n;
        fscanf(file, "%d\n", &n);
        if (n > 0) {
            mat = create_empty(n);
            if (mat) {
                for (int i=0; i<n; i++) {
                    for (int j=0; j<n + 1; j++) {
                        if (sizeof(M_TYPE) == sizeof(float)) {
                            fscanf(file, "%f", mat->data[i] + j);
                        }
                        if (sizeof(M_TYPE) == sizeof(double)) {
                            fscanf(file, "%lf", mat->data[i] + j);
                        }
                    }
                }
            }
        }
    } else {
        printf("File not found!\n");
    }
    if (file) fclose(file);

    return mat;
}

void delete_matrix(matrix_t * matrix) {
    if (matrix) {
        for (int i=0; i<matrix->n; i++)
            free(matrix->data[i]);
        free(matrix->data);
        free(matrix);
    }
}

void spaw_rows(int i, int j, matrix_t *matrix) {
    M_TYPE *tmp = matrix->data[i];
    matrix->data[i] = matrix->data[j];
    matrix->data[j] = tmp;
}

int get_row_max_col(matrix_t *mat, int idx) {
    int max_idx = -1;
    M_TYPE max_num = 0;
    for (int it=idx; it<mat->n; it++) {
        if (ABS(mat->data[it][idx]) > ABS(max_num)) {
            max_num = mat->data[it][idx];
            max_idx = it;
        }
    }
    return max_idx;
}

// returns presicion in form of decimal places
int calculate_diff(int n, M_TYPE *orig_sols, M_TYPE *sols) {
    M_TYPE diff;
    M_TYPE max_diff = 0;
    for (int i=0; i<n; i++) {
        diff = ABS(sols[i] - orig_sols[i]);
        max_diff = diff > max_diff ? diff : max_diff;
    }
    int p = 0;
    while (max_diff < 1 && max_diff) {
        p++;
        max_diff *= 10;
    }
    return p;
}

// returns solutions
M_TYPE *Gauss_solve(matrix_t *matt) {
    int m_idx;
    int broke = 0;
    M_TYPE *sols;

    if (matt) {
        sols = (M_TYPE *) calloc(matt->n, sizeof(M_TYPE));

        // direct move
        for (int it=0; it<matt->n; it++) {
            // swap row with max column element
            m_idx = get_row_max_col(matt, it);
            if (m_idx == -1) {
                printf("Infinite solutions case! Over\n");
                broke = 1;
                break;
            }
            if (m_idx != it) spaw_rows(it, m_idx, matt);
            // substract from each
            for (int s_it = it + 1; s_it<matt->n; s_it++) {
                M_TYPE sub = matt->data[s_it][it] / matt->data[it][it];
                for (int k=it; k<matt->n + 1; k++) {
                    matt->data[s_it][k] -= sub * matt->data[it][k];
                }
            }
        }

        // reverse move
        for (int i=matt->n - 1; (i>=0) && (! broke); i--) {
            M_TYPE sum = 0;
            for (int j=i + 1; j<matt->n; j++) {
                sum += matt->data[i][j] * sols[j];
            }
            sols[i] = (matt->data[i][matt->n] - sum) / matt->data[i][i];
        }
    }

    return sols;
}

void step(int n, FILE* file) {
    clock_t start, end;
    double dur;
    M_TYPE *orig = (M_TYPE *) calloc(n, sizeof(M_TYPE));
    matrix_t *matrix = generate_random_sols(n, orig);
    start = clock();
    M_TYPE *sols = Gauss_solve(matrix);
    // save time in ms
    end = clock();
    dur = (end - start) * 1000 * 1000 / CLOCKS_PER_SEC;
    float pr = calculate_diff(n, orig, sols);
    fprintf(file, "%d %lf %f\n", n, dur, pr);
    free(orig);
    free(sols);
    delete_matrix(matrix);
}

int main(int argc, char * argv[]) {
    srand(time(NULL));
    matrix_t *matrix = 0;
    // explore complexity and precision
    // from 1 to 9
    FILE *file = fopen(argv[1], "w");
    for (int n=1; n<10; n++) {
        step(n, file);
    }
    // from 10 to 100 with step 10
    for (int n=10; n<100; n += 10) {
        step(n, file);
    }
    // from 100 to 10'000 with step 100
    for (int n=100; n<=5000; n += 100) {
        step(n, file);
    }

    return 0;
}
