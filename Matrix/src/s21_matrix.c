#include "s21_matrix.h"

#include <math.h>
#include <stdlib.h>

int s21_create_matrix(int rows, int columns, matrix_t *result) {
    int error = 0;
    result->rows = rows;
    result->columns = columns;
    if (rows >= 1 && columns >= 1) {
        result->matrix = (double**)calloc(rows, sizeof(double*));
            for (int i = 0; i < rows; i++) {
                result->matrix[i] = (double*)calloc(columns, sizeof(double));
                for (int j = 0; j < columns; j++)
                    result->matrix[i][j] = 0.0;
            }
    } else {
        error = 1;
        result->matrix = NULL;
    }
    return error;
}

void s21_remove_matrix(matrix_t *A) {
    if (A->matrix != NULL) {
        for (int i = 0; i < A->rows; i++)
            free(A->matrix[i]);
        free(A->matrix);
        A->matrix = NULL;
    }
    if (A->rows)
        A->rows = 0;
    if (A->columns)
        A->columns = 0;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
    int error = SUCCESS;
    if ((A->rows != B->rows) || (A->columns != B->columns))
        error = FAILURE;
    for (int i = 0; i < A->rows && error == SUCCESS; i++)
        for (int j = 0; j < B->columns && error == SUCCESS; j++)
            if (fabs(A->matrix[i][j] - B->matrix[i][j]) > EPS)
                error = FAILURE;
    return error;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
    int error = 0;
    if (A->rows < 1 || B->rows < 1 || A->columns < 1 || B->columns < 1
    || A->matrix == NULL || B->matrix == NULL) {
        error = 1;
    } else if (A->rows != B->rows || A->columns != B->columns) {
        error = 2;
    } else {
        s21_create_matrix(A->rows, A->columns, result);
        for (int i = 0; i < A->rows; i++)
            for (int j = 0; j < A->columns; j++)
                result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
    }
    return error;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
    int error = 0;
    if (A->rows < 1 || B->rows < 1 || A->columns < 1 || B->columns < 1
    || A->matrix == NULL || B->matrix == NULL) {
        error = 1;
    } else if (A->rows != B->rows || A->columns != B->columns) {
        error = 2;
    } else {
        s21_create_matrix(A->rows, A->columns, result);
        for (int i = 0; i < A->rows; i++)
            for (int j = 0; j < A->columns; j++)
                result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
    }
    return error;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
    int error = 0;
    if (A->columns < 1 || A->rows < 1 || A->matrix == NULL) {
        error = 1;
    } else {
        s21_create_matrix(A->rows, A->columns, result);
        for (int i = 0; i < A->rows; i++)
            for (int j = 0; j < A->columns; j++)
                result->matrix[i][j] = A->matrix[i][j] * number;
    }
    return error;
}


int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
    int error = 0;
    if (A->rows < 1 || B->rows < 1 || A->columns < 1 || B->columns < 1
    || A->matrix == NULL || B->matrix == NULL) {
        error = 1;
    } else if (B->rows != A->columns) {
        error = 2;
    } else {
        s21_create_matrix(A->rows, B->columns, result);
        for (int i = 0; i < A->rows; i++)
            for (int j = 0; j < B->columns; j++) {
                int summa = 0;
                for (int k = 0; k < A->columns; k++)
                    summa += A->matrix[i][k] * B->matrix[k][j];
                result->matrix[i][j] = summa;
            }
    }
    return error;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
    int error = 0;
    if (A->rows < 1 || A->columns < 1 || A->matrix == NULL) {
        error = 1;
    } else {
        s21_create_matrix(A->columns, A->rows, result);
        for (int i = 0; i < A->columns; i++)
            for (int j = 0; j < A->rows; j++)
                result->matrix[i][j] = A->matrix[j][i];
    }
    return error;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
    int error = 0;
    if (A != NULL && A->matrix != NULL && A->rows > 0 && A->columns > 0) {
        if (A->columns == A->rows) {
            s21_create_matrix(A->rows, A->columns, result);
            for (int i = 0; i < A->rows; i++) {
                for (int j = 0; j < A->columns; j++) {
                    double determinant = 0;
                    matrix_t tempMatrix;
                    s21_create_matrix(A->rows - 1, A->columns - 1, &tempMatrix);
                    get_temp_matrix(A, &tempMatrix, i, j);
                    s21_determinant(&tempMatrix, &determinant);
                    result->matrix[i][j] = pow(-1, i + j) * determinant;
                    s21_remove_matrix(&tempMatrix);
                }
            }
        } else {
            error = 2;
        }
    } else {
        error = 1;
    }
    return error;
}

int s21_determinant(matrix_t *A, double *result) {
    int error = 0;
    if (A->rows < 1 || A->columns < 1 || !A->matrix) {
        error = 1;
    } else if (A->rows != A->columns) {
        error = 2;
    } else {
        if (A->rows == 1) {
            *result = A->matrix[0][0];
        } else if (A->rows == 2) {
            *result = A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
        } else {
            double determinant = 0;
            for (int j = 0; j < A->columns; j++) {
                matrix_t tempMatrix;
                s21_create_matrix(A->rows - 1, A->rows - 1, &tempMatrix);
                get_temp_matrix(A, &tempMatrix, 0, j);
                s21_determinant(&tempMatrix, result);
                determinant += pow(-1, j) * A->matrix[0][j] * (*result);
                s21_remove_matrix(&tempMatrix);
            }
            *result = determinant;
        }
    }
    return error;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
    int error = 0;
    if (A->rows < 1 || A->columns < 1 || A->matrix == NULL) {
        error = 1;
    } else if (A->rows != A->columns) {
        error = 2;
    } else {
        double determinant = 0.0;
        s21_determinant(A, &determinant);
        if (fabs(determinant) > EPS) {
            matrix_t tempMatrix1, tempMatrix2;
            s21_calc_complements(A, &tempMatrix1);
            s21_transpose(&tempMatrix1, &tempMatrix2);
            s21_mult_number(&tempMatrix2, 1.0 / determinant, result);
            s21_remove_matrix(&tempMatrix1);
            s21_remove_matrix(&tempMatrix2);
        } else {
            error = 2;
        }
    }
    return error;
}

void get_temp_matrix(matrix_t *A, matrix_t *tempMatrix, int row, int column) {
    int tempMatrix_i = 0, tempMatrix_j = 0;
    for (int i = 0; i < A->rows; i++)
        for (int j = 0; j < A->columns; j++)
            if (i != row && j != column) {
                tempMatrix->matrix[tempMatrix_i][tempMatrix_j] = A->matrix[i][j];
                if (tempMatrix_j + 1 != tempMatrix->columns) {
                    tempMatrix_j++;
                } else {
                    tempMatrix_i++;
                    tempMatrix_j = 0;
                }
            }
}
