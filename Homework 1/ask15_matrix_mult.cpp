#include "matrix_mult.h"

int cartesian2flat_row_major(int i, int j, int n_col);
int cartesian2flat_col_major(int i, int j, int n_row);
void matrix_vec_mult_row_major(double** M, double* v, double* result, int n_row, 
int n_col);
void matrix_vec_mult_col_major(double** M, double* v, double* result, int n_row, 
int n_col);
void matrix_vec_mult_sparse(SparseMatrix M, double* v, double* result, int n);

int cartesian2flat_row_major(int i, int j, int n_col) {
    return i*n_col+j;
}

int cartesian2flat_col_major(int i, int j, int n_row) {
    return i*n_row+j;
}

void matrix_vec_mult_row_major(double** M, double* v, double* result, int n_row, int n_col) {
    for (int i = 0; i < n_row; i++) {
        double total = 0;
        for (int j = 0; j < n_col; j++) {
            total += M[i][j]*v[j]; //inner product of row vectors of the matrix by the column vector
        }
        result[i] = total;
    }
}

void matrix_vec_mult_col_major(double** M, double* v, double* result, int n_row, int n_col) {
    for (int j = 0; j < n_col; j++) {
        for (int i = 0; i < n_row; i++) {
            result[i] += M[i][j]*v[j]; //a linear combination of the columns of the matrix
        }
    }
}

void matrix_vec_mult_sparse(SparseMatrix M, double* v, double* result, int n) { //n  is the size of col_ptr array
    for (int i = 1; i <= n; i++) {
        for (int j = M.col_ptr[i-1]; j < M.col_ptr[i]; j++) {
            result[M.row_ind[j]] += M.val[j]*v[i-1]; //getting the value from the val matrix corresponding to the right row index and right column before multiplying it with the value in the vector
        }
    }
}