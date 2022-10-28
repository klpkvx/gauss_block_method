#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include "time.h"
#include "init.h"
#include "solution.h"
#include "matrix.h"
#include "errors.h"

int main(int argc, char **argv)
{
  FILE *file = nullptr;
  int task = 13;
  int n = 0, m = 0, r = 0, formulae = 0;
  int *index = nullptr;
  double *matrix = nullptr;         // Исходная матрица
  double *inverse_matrix = nullptr; // Обратная матрица
  double *result_matrix = nullptr;  // для расчета невязки
  double *block1 = nullptr, *block2 = nullptr, *block3 = nullptr, *block4 = nullptr;
  double *unit_matrix = nullptr; // Дополнительная матрица для расчета невязки
  double t1 = 0, t2 = 0;
  double r1 = 0, r2 = 0;
  double matrix_norm = 0;
  char *filename = nullptr;
  int res = 0, result_calc = 0;

  res = check_args(argc, argv, &filename, n, m, r, formulae);
  if (res != (int)ERRORS::SUCCESS)
    return res;

  if (argc == 6)
  {
    file = fopen(filename, "rt");
    if (!file)
    {
      fprintf(stderr, "Cannot open file!\n");
      return (int)ERRORS::CANNOT_OPEN;
    }
  }

  matrix = new double[n * n];
  if (!matrix)
  {
    fprintf(stderr, "Cannot allocate memory!\n");
    return (int)ERRORS::CANNOT_ALLOCATE;
  }

  res = init_matrix(&matrix, formulae, n, file, argc);
  if (res != (int)ERRORS::SUCCESS)
  {
    delete[] matrix;
    if (formulae == 0)
      fclose(file);
    fprintf(stdout, "CODE:res %d!\n", res);
    return -1;
  }
  matrix_norm = norma(matrix, n);

  fprintf(stdout, "matrix\n");
  print_matrix(r, n, matrix);
  fprintf(stdout, "\n");

  inverse_matrix = new double[n * n];
  if (!inverse_matrix)
  {
    delete[] matrix;
    fprintf(stderr, "Cannot allocate memory!\n");
    return (int)ERRORS::CANNOT_ALLOCATE;
  }

  init_inverse_matrix(&inverse_matrix, n);

  index = new int[m];
  if (!index)
  {
    delete[] matrix;
    delete[] inverse_matrix;
    return -1;
  }
  for (int i = 0; i < m; i++)
    index[i] = i;

  int *index_block = new int[n];
  for (int i = 0; i < n; i++)
    index_block[i] = i;

  block1 = new double[m * m];
  block2 = new double[m * m];
  block3 = new double[m * m];
  block4 = new double[m * m];

  for (int i = 0; i < m * m; i++)
  {
    block1[i] = 0;
    block2[i] = 0;
    block3[i] = 0;
    block4[i] = 0;
  }

  // Поиск обратной матрицы
  t1 = clock();
  result_calc = gauss_method(matrix, inverse_matrix, index, index_block, block1, block2, block3, block4, n, m, matrix_norm);
  t1 = (clock() - t2) / CLOCKS_PER_SEC;

  if (formulae == 0)
    rewind(file);

  res = init_matrix(&matrix, formulae, n, file, argc);
  if (res != (int)ERRORS::SUCCESS)
  {
    delete[] matrix;
    delete[] inverse_matrix;
    if (formulae == 0)
      fclose(file);
    fprintf(stdout, "CODE:res %d!\n", res);
    return -1;
  }
  if (formulae == 0)
    fclose(file);

  unit_matrix = new double[n * n];
  if (!unit_matrix)
  {
    fprintf(stderr, "Cannot allocate!\n");
    delete[] matrix;
    delete[] inverse_matrix;
    if (formulae == 0)
      fclose(file);
    return (int)ERRORS::CANNOT_ALLOCATE;
  }
  else
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        if (i == j)
          unit_matrix[i * n + j] = 1.;
        else
          unit_matrix[i * n + j] = 0.;
      }
    }

  result_matrix = new double[n * n]; // Для расчета невязки
  if (!result_matrix)
  {
    fprintf(stderr, "Cannot allocate!\n");
    delete[] matrix;
    delete[] inverse_matrix;
    delete[] unit_matrix;
    if (formulae == 0)
      fclose(file);
    return (int)ERRORS::CANNOT_ALLOCATE;
  }

  for (int i = 0; i < n * n; i++)
    result_matrix[i] = 0;

  t2 = clock(); // Вычисление невязки
  residual(r1, r2, matrix, inverse_matrix, result_matrix, unit_matrix, n, m, result_calc, r);
  t2 = (clock() - t2) / CLOCKS_PER_SEC;

  printf("\n");
  printf(
      "%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n",
      argv[0], task, r1, r2, t1, t2, formulae, n, m);
  delete[] matrix;
  delete[] inverse_matrix;
  delete[] result_matrix;
  delete[] unit_matrix;
  delete[] index;
  delete[] block1;
  delete[] block2;
  delete[] block3;
  delete[] block4;
  delete[] index_block;
  return 0;
}
