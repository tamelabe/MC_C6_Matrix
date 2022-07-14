#ifndef SRC_S21_MATRIX_H_
#define SRC_S21_MATRIX_H_
#define OK 0
#define INCRT_MTRX 1
#define CALC_ERROR 2
#define SUCCESS 1
#define FAILURE 0

typedef struct matrix_struct {
    double** matrix;
    int rows;
    int columns;
} matrix_t;


#endif  // SRC_S21_MATRIX_H_

