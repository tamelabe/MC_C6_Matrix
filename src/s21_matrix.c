#include "s21_matrix.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// BASIC FUNCTIONS
int s21_create_matrix(int rows, int columns, matrix_t *result);     // done
void s21_remove_matrix(matrix_t *A);                                // done
int s21_eq_matrix(matrix_t *A, matrix_t *B);                        // done
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);     // done
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);     // done
int s21_mult_number(matrix_t *A, double number, matrix_t *result);  // done
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);    // done
int s21_transpose(matrix_t *A, matrix_t *result);                   // done
int s21_calc_complements(matrix_t *A, matrix_t *result);
int s21_determinant(matrix_t *A, double *result);
int s21_inverse_matrix(matrix_t *A, matrix_t *result);

// ADDITIONAL FUNCTIONS
int s21_check_mex(matrix_t *A);
int s21_check_mrc(matrix_t *A);
int s21_check_msz(matrix_t *A, matrix_t *B);
void s21_rm_rc(matrix_t *BF, int nj, matrix_t *AF);

int main() {
  int rows = 3, columns = 3;
  double det = 0;
  int cnt_test = 0;
  matrix_t A;  //, B, result;
  s21_create_matrix(3, 3, &A);
  // s21_create_matrix(2, 3, &B);
  printf("Fill A matrix: ");
  for (int i = 0; i < A.rows; i++) {
    for (int j = 0; j < A.columns; j++) {
      A.matrix[i][j] = cnt_test;
      cnt_test++;
      // scanf("%lf", &A.matrix[i][j]);
    }
  }
  /*
  printf("Fill B matrix: ");
  for (int i = 0; i < B.rows; i++) {
    for (int j = 0; j < B.columns; j++) {
      scanf("%lf", &B.matrix[i][j]);
    }
  }
  */
  // s21_mult_matrix(&A, &B, &result);
  // s21_remove_matrix(&A);
  // s21_remove_matrix(&B);
  // printf("Mult is: \n");
  // for (int i = 0; i < columns; i++) {
  //   for (int j = 0; j < rows; j++) {
  //     printf("%-2.0lf", result.matrix[i][j]);
  //     if (j == columns - 1 && i < rows - 1) {
  //       putchar('\n');
  //     } else if (j != columns - 1) {
  //       putchar(' ');
  //     }
  //   }
  // }
  s21_determinant(&A, &det);
  s21_remove_matrix(&A);
  printf("Determinant is: %.0lf\n", det);
  return 0;
}

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int res_code = INCRT_MTRX;
  result->matrix = NULL;
  if (rows > 0 && columns > 0) {
    result->matrix = calloc(rows, sizeof(double *));
    if (result->matrix != NULL) {
      res_code = OK;
      for (int i = 0; i < rows; i++) {
        result->matrix[i] = calloc(columns, sizeof(double));
      }
    }
    result->rows = rows;
    result->columns = columns;
  }
  return res_code;
}

void s21_remove_matrix(matrix_t *A) {
  if (A->matrix != NULL) {
    for (int i = 0; i < A->rows; i++) {
      free(A->matrix[i]);
    }
    free(A->matrix);
    A->rows = 0;
    A->columns = 0;
    A = NULL;
  }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int res_code = FAILURE;
  if (s21_check_mex(A) && s21_check_mex(B) && s21_check_msz(A, B)) {
    res_code = SUCCESS;
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; i++) {
        if (fabs(A->matrix[i][j] - B->matrix[i][j] > 1e-7)) {
          res_code--;
        }
      }
    }
  }
  if (res_code < 1) res_code = FAILURE;
  return res_code;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res_code = OK;
  if (!s21_check_mex(A) || !s21_check_mex(B)) {
    res_code = INCRT_MTRX;
  } else {
    if (!s21_check_msz(A, B)) {
      res_code = CALC_ERROR;
    } else if (!(res_code = s21_create_matrix(A->rows, A->columns, result))) {
      for (int i = 0; i < A->rows; i++) {
        for (int j; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
        }
      }
    }
  }
  return res_code;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res_code = OK;
  if (!s21_check_mex(A) || !s21_check_mex(B)) {
    res_code = INCRT_MTRX;
  } else {
    if (!s21_check_msz(A, B)) {
      res_code = CALC_ERROR;
    } else if (!(res_code = s21_create_matrix(A->rows, A->columns, result))) {
      for (int i = 0; i < A->rows; i++) {
        for (int j; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
        }
      }
    }
  }
  return res_code;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int res_code = OK;
  if (!s21_check_mex(A)) {
    res_code = INCRT_MTRX;
  } else if (!(res_code = s21_create_matrix(A->rows, A->columns, result))) {
    for (int i = 0; i < A->rows; i++) {
      for (int j; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] * number;
      }
    }
  }
  return res_code;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res_code = OK;
  if (!s21_check_mex(A) || !s21_check_mex(B)) {
    res_code = INCRT_MTRX;
  } else if (A->columns == B->rows) {
    if (!(res_code = s21_create_matrix(A->rows, B->columns, result))) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < B->columns; j++) {
          for (int n = 0; n < B->rows; n++) {
            result->matrix[i][j] += A->matrix[i][n] * B->matrix[n][j];
          }
        }
      }
    }
  } else {
    res_code == CALC_ERROR;
  }
  return res_code;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int res_code = INCRT_MTRX;
  if (!s21_check_mex(A) || !s21_check_mrc(A)) {
    res_code = INCRT_MTRX;
  } else if (!(res_code = s21_create_matrix(A->columns, A->rows, result))) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[j][i];
      }
    }
  }
  return res_code;
}

int s21_determinant(matrix_t *A, double *result) {
  int res_code = OK;
  double degree = 1.0;
  if (!s21_check_mex(A) || !s21_check_mrc(A)) {
    res_code = INCRT_MTRX;
  } else if (A->rows != A->columns) {
    res_code = CALC_ERROR;
  } else {
    if (A->rows == 1) {
      *result = A->matrix[0][0];
    } else if (A->rows == 2) {
      *result == (A->matrix[0][0] * A->matrix[1][1]) -
                     (A->matrix[0][1] * A->matrix[1][0]);
    } else {
      matrix_t Rmatrix;
      if (!(res_code =
                s21_create_matrix(A->columns - 1, A->rows - 1, &Rmatrix))) {
        for (int j = 0; j < A->rows; j++) {
          s21_rm_rc(A, j, &Rmatrix);
          *result +=
              (degree * A->matrix[0][j] + s21_determinant(&Rmatrix, result));
          degree = -degree;
          printf("result = %.0lf\n", *result);
        }
        for (int i = 0; i < Rmatrix.columns; i++) {
          for (int j = 0; j < Rmatrix.rows; j++) {
            printf("%-2.0lf", Rmatrix.matrix[i][j]);
            if (j == Rmatrix.rows - 1 && i < Rmatrix.columns - 1) {
              putchar('\n');
            } else if (j != Rmatrix.columns - 1) {
              putchar(' ');
            }
          }
        }
        s21_remove_matrix(&Rmatrix);
      }
    }
  }
  return res_code;
}

void s21_rm_rc(matrix_t *BF, int nj, matrix_t *AF) {
  int irow = 0, icol = 0;
  for (int i = 0; i < BF->rows - 1; i++) {
    if (!i) {
      irow = 1;
    }
    for (int j = 0; j < BF->rows - 1; j++) {
      if (j == nj) {
        icol = 1;
      }
      AF->matrix[i][j] = BF->matrix[i + irow][j + icol];
    }
  }
}

int s21_check_mex(matrix_t *A) {
  int flag = FAILURE;
  if (A != NULL && A->matrix != NULL) {
    flag = SUCCESS;
  }
  return flag;
}

int s21_check_mrc(matrix_t *A) {
  int flag = FAILURE;
  if (A->rows > 0 && A->columns > 0) {
    flag = SUCCESS;
  }
  return flag;
}

int s21_check_msz(matrix_t *A, matrix_t *B) {
  int flag = FAILURE;
  if (A->rows == B->rows && A->columns == B->columns) {
    flag = SUCCESS;
  }
  return flag;
}
