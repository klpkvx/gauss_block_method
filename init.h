#ifndef INIT_H
#define INIT_H
#include "stdio.h"

double f(int formulae, int n, int i, int j);
int max(int i, int j);

int init_matrix(double **matrix, int formulae, int n, FILE *input, int argc);
void init_inverse_matrix(double **attached_matrix, int n);

#endif // INIT_H
