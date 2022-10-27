#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "init.h"
#include "errors.h"

void init_inverse_matrix(double **attached_matrix, int n)
{
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
    {
      if (i == j)
        (*attached_matrix)[i * n + j] = 1;
      else
        (*attached_matrix)[i * n + j] = 0;
    }
}

int init_matrix(double **matrix, int formulae, int n, FILE *input, int argc)
{

  if (argc == 5)
  {
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        (*matrix)[i * n + j] = f(formulae, n, i, j);
  }
  else
  {
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        if (fscanf(input, "%lf", &(*matrix)[i * n + j]) != 1)
        {
          if (!feof(input))
          {
            fprintf(stdout,
                    "Cannot read element a[%d,%d] in file \n", i,
                    j);
          }
          else
            fprintf(stdout, "Not enough elements to read !\n");
          return -4;
        }
  }

  return (int)ERRORS::SUCCESS;
}

int max(int i, int j)
{
  return i > j ? i : j;
}

double
f(int formulae, int n, int i_orig, int j_orig)
{
  int i = i_orig + 1;
  int j = j_orig + 1;
  switch (formulae)
  {
  case 1:
    return n - max(i, j) + 1;
  case 2:
    return max(i, j);
  case 3:
    return fabs(i - j);
  case 4:
    return 1. / (i + j - 1);
  }
  return 0;
}
