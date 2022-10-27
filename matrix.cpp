#include "stdlib.h"
#include "stdio.h"
#include "solution.h"
#include "matrix.h"
#include "errors.h"

int check_args(int &argc, char **argv, char **filename, int &n, int &m, int &r,
							 int &formulae)
{
	if (argc == 5 || argc == 6)
	{
		if (argc == 5)
		{
			if (!((sscanf(argv[1], "%d", &n) == 1) && (n > 0) && (sscanf(argv[2], "%d", &m) == 1) && (m > 0) && (sscanf(argv[3], "%d", &r) == 1) && (r >= 0) && (sscanf(argv[4], "%d", &formulae) == 1) && (formulae > 0 && formulae <= 4)))
			{
				fprintf(stdout, "Usage: %s n m r s\n", argv[0]);
				return (int)ERRORS::ARGS;
			}
		}
		else
		{
			if (argc == 6)
			{
				if (!((sscanf(argv[1], "%d", &n) == 1) && (n > 0) && (sscanf(argv[2], "%d", &m) == 1) && (m > 0) && (sscanf(argv[3], "%d", &r) == 1) && (r >= 0) && (sscanf(argv[4], "%d", &formulae) == 1) && (formulae == 0)))
				{
					fprintf(stdout, "Usage: %s n m r 0 filename\n", argv[0]);
					return (int)ERRORS::ARGS;
				}
				*filename = argv[5];
				if (!filename)
					return (int)ERRORS::ARGS;
			}
		}
	}
	else
	{
		fprintf(stdout, "Usage: %s n m r s filename\n", argv[0]);
		return (int)ERRORS::ARGS;
	}

	if (m > n)
	{
		fprintf(stderr, "m > n!\n");
		return (int)ERRORS::ARGS;
	}

	return (int)ERRORS::SUCCESS;
}

void print_block(double *block, int row, int col)
{

	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
			printf("%10.3e ", block[i * col + j]);
		printf("\n");
	}
}

void print_matrix(int r, int n, double *array)
{
	r < n ? r : r = n;
	for (int i = 0; i < r; i++)
	{
		{
			for (int j = 0; j < r; j++)
				printf(" %10.3e", array[i * n + j]);
		}
		printf("\n");
	}
}
