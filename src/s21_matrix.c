#include "s21_matrix.h"

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
      for (int j = 0; j < A->columns; j++) {
        if (fabs(A->matrix[i][j] - B->matrix[i][j]) > 1e-7) {
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
  } else if (!s21_check_msz(A, B)) {
      res_code = CALC_ERROR;
  } else if (!(res_code = s21_create_matrix(A->rows, A->columns, result))) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
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
        for (int j = 0; j < A->columns; j++) {
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
      for (int j = 0; j < A->columns; j++) {
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
    res_code = CALC_ERROR;
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
        result->matrix[j][i] = A->matrix[i][j];
      }
    }
  }
  return res_code;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int res_code = OK;
  if (!s21_check_mex(A) || !s21_check_mrc(A)) {
    res_code = INCRT_MTRX;
  } else if (A->rows != A->columns) {
    res_code = CALC_ERROR;
  } else {
    if (!(res_code = s21_create_matrix(A->rows, A->columns, result))) {
      res_code = s21_convert_to_minor(A, result);
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] *= pow(-1, i + j);
        }
      }
    }
  }
  return res_code;
}

int s21_determinant(matrix_t *A, double *result) {
  int res_code = OK;
  double degree = 1.0, nres = 0;
  if (!s21_check_mex(A) || !s21_check_mrc(A)) {
    res_code = INCRT_MTRX;
  } else if (A->rows != A->columns) {
    res_code = CALC_ERROR;
  } else {
    if (A->rows == 1) {
      *result = A->matrix[0][0];
    } else if (A->rows == 2) {
      *result = (A->matrix[0][0] * A->matrix[1][1]) -
                (A->matrix[0][1] * A->matrix[1][0]);
    } else {
      matrix_t Rmatrix;
      if (!(res_code =
                s21_create_matrix(A->columns - 1, A->rows - 1, &Rmatrix))) {
        for (int j = 0; j < A->rows; j++) {
          s21_rm_rc(A, 0, j, &Rmatrix);
          s21_determinant(&Rmatrix, &nres);
          *result += degree * A->matrix[0][j] * nres;
          degree = -degree;
          nres = 0;
        }
        s21_remove_matrix(&Rmatrix);
      }
    }
  }
  return res_code;
}

int s21_convert_to_minor(matrix_t *A, matrix_t *MM) {
  int res_code = OK;
  matrix_t buf;
  if (!(res_code = s21_create_matrix(A->rows - 1, A->columns - 1, &buf))) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        double det = 0.0;
        s21_rm_rc(A, i, j, &buf);
        s21_determinant(&buf, &det);
        MM->matrix[i][j] = det;
      }
    }
    s21_remove_matrix(&buf);
  }
  return res_code;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int res_code = OK;
  if (!s21_check_mex(A) || !s21_check_mrc(A)) {
    res_code = INCRT_MTRX;
  } else if (A->rows != A->columns) {
    res_code = CALC_ERROR;
  } else {
    double det = 0;
    s21_determinant(A, &det);
    if (fabs(det) > 1e-7) {
      matrix_t tempM1, tempM2;
      s21_calc_complements(A, &tempM1);
      s21_transpose(&tempM1, &tempM2);
      s21_remove_matrix(&tempM1);
      s21_mult_number(&tempM2, 1.0 / det, result);
      s21_remove_matrix(&tempM2);
    } else {
      res_code = CALC_ERROR;
    }
  }
  return res_code;
}

void s21_rm_rc(matrix_t *BF, int ni, int nj, matrix_t *AF) {
  int irow = 0, icol = 0;
  for (int i = 0; i < AF->rows; i++) {
    if (i == ni) {
      irow = 1;
    }
    for (int j = 0; j < AF->columns; j++) {
      if (j == nj) {
        icol = 1;
      }
      AF->matrix[i][j] = BF->matrix[i + irow][j + icol];
    }
    icol = 0;
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
