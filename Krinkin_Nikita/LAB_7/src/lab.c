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


int main(int argc, char *argv[]) {
    srand(time(NULL));

    int n;
    matrix_t *matrix = 0;
    if (argc == 1) {
        printf("Generating random table\n");
        printf("intput table size: ");
        scanf("%d", &n);
        if (n <= 0) {
            printf("Incorrect size!\n");
            goto end;
        }
        matrix = generate_random(n);
    }

    if (argc == 2) {
        printf("Getting table from file\n");
        matrix = from_file(argv[1]);
    }

    if (matrix) {
        printf("Your table: \n");
        print_table(matrix);
    }

    if (matrix) {
        M_TYPE *sols = Gauss_solve(matrix);
        if (sols) {
            for (int i=0; i<matrix->n; i++) {
                if (sizeof(M_TYPE) == sizeof(float)) printf("X%d = %f\n", i, sols[i]);
                if (sizeof(M_TYPE) == sizeof(double)) printf("X%d = %lf\n", i, sols[i]);
            }
        }
    }

end:
    if (matrix) delete_matrix(matrix);
    return 0;
}
