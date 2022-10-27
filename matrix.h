#ifndef MATRIX_H
#define MATRIX_H

int check_args(int &argc, char **argv, char **filename, int &n, int &m, int &r, int &formulae);
void print_block(double *block, int row, int col);
void print_matrix_blocks(double **matrix, double **block, int n, int m);
void print_matrix(int r, int n, double *array);
#endif // MATRIX_H