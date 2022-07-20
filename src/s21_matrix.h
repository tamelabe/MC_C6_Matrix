#ifndef SRC_S21_MATRIX_H_
#define SRC_S21_MATRIX_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define OK 0
#define INCRT_MTRX 1
#define CALC_ERROR 2
#define SUCCESS 1
#define FAILURE 0

typedef struct matrix_struct {
  double **matrix;
  int rows;
  int columns;
} matrix_t;

// BASIC FUNCTIONS
int s21_create_matrix(int rows, int columns, matrix_t *result);
void s21_remove_matrix(matrix_t *A);
int s21_eq_matrix(matrix_t *A, matrix_t *B);
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_mult_number(matrix_t *A, double number, matrix_t *result);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_transpose(matrix_t *A, matrix_t *result);
int s21_calc_complements(matrix_t *A, matrix_t *result);
int s21_determinant(matrix_t *A, double *result);
int s21_inverse_matrix(matrix_t *A, matrix_t *result);
// ADDITIONAL FUNCTIONS
int s21_check_mex(matrix_t *A);
int s21_check_mrc(matrix_t *A);
int s21_check_msz(matrix_t *A, matrix_t *B);
void s21_rm_rc(matrix_t *BF, int ni, int nj, matrix_t *AF);
int s21_convert_to_minor(matrix_t *A, matrix_t *MM);

#endif  // SRC_S21_MATRIX_H_
