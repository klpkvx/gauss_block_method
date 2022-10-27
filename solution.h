#ifndef SOLUTION_H
#define SOLUTION_H

void get_block(double *matrix, double *block, int i_block, int j_block, int n, int m);
void set_block(double *matrix, double *block, int i_block, int j_block, int n, int m);
void set_unit_zero_block(double *matrix, double *block, int i_block, int j_block, int n, int m);
double norma(double *matrix, int n);

void residual(double &r1, double &r2, double *matrix, double *inverse_matrix,
							double *block, double *sum_array, int n, int m, int res, int r);
double find_residual(double *matrix, double *inversed_matrix, double *block, double *sum_array, int n, int m);

int gauss_classic_row(double *matrix, double *inverse_matrix, int *index, int n, double matrix_norm, int row_ind);
int gauss_method(double *matrix, double *inversed_matrix, int *index, int *index_block, double *block1, double *block2, double *block3, double *block4, int n, int m, double matrix_norm);
void mult(double *pa, double *pb, double *pc, int m1, int m2, int m3, int m);
void matrix_minus(double *a, double *b, int i_block, int j_block, int m, int k, int l);
#endif // SOLUTION_H