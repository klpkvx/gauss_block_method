#include "stdio.h"
#include "solution.h"
#include "memory.h"
#include "math.h"
#include "errors.h"
#include "matrix.h"

#define EPS 1e-16
#define CANNOT_SOLVE -1
#define IRREVERSIBLE -2
#define SUCCESS 0

void get_block(double *matrix, double *block, int i_block, int j_block, int n, int m)
{
	int k = n / m;
	int l = n % m;
	int h = (i_block < k ? m : l);
	int v = (j_block < k ? m : l);
	int tmp;
	for (int i = 0; i < m; i++)
	{
		tmp = i * m;
		for (int j = 0; j < m; j++)
			if (i < h && j < v)
				block[tmp + j] = matrix[n * (i_block * m + i) + j_block * m + j];
			else
				block[tmp + j] = 0;
	}
}

void set_block(double *matrix, double *block, int i_block, int j_block, int n, int m)
{
	int k = n / m;
	int l = n % m;
	int h = (i_block < k ? m : l);
	int v = (j_block < k ? m : l);
	int tmp;
	for (int i = 0; i < h; i++)
	{
		tmp = i * m;
		for (int j = 0; j < v; j++)
			matrix[n * (i_block * m + i) + j_block * m + j] = block[tmp + j];
	}
}

double norma(double *matrix, int n)
{
	double norm = 0, tmp_sum = 0;
	for (int i = 0; i < n; i++)
	{
		tmp_sum = 0;
		for (int j = 0; j < n; j++)
			tmp_sum += fabs(matrix[j * n + i]);
		if (tmp_sum > norm)
			norm = tmp_sum;
	}
	return norm;
}

void mult(double *a, double *b, double *res, int m1, int m2, int m3, int m)
{
	// a m1 x m2
	// b m2 x m3
	int t = 0, q = 0, r = 0;
	int v = m1, h = m3, ah = m2;
	int v3 = v % 3, h3 = h % 3;
	double s00 = 0, s01 = 0, s02 = 0;
	double s10 = 0, s11 = 0, s12 = 0;
	double s20 = 0, s21 = 0, s22 = 0;

	for (r = 0; r < v; r++)
		for (t = 0; t < h; t++)
			res[r * m + t] = 0;

	for (r = 0; r < v3; r++)
	{
		for (t = 0; t < h3; t++)
		{
			double sum;
			sum = 0;
			for (q = 0; q < ah; q++)
				sum += a[r * m + q] * b[q * m + t];
			res[r * m + t] += sum;
		}
		for (; t < h; t += 3)
		{
			s00 = 0;
			s01 = 0;
			s02 = 0;
			for (q = 0; q < ah; q++)
			{
				s00 += a[r * m + q] * b[q * m + t];
				s01 += a[r * m + q] * b[q * m + t + 1];
				s02 += a[r * m + q] * b[q * m + t + 2];
			}
			res[r * m + t] += s00;
			res[r * m + t + 1] += s01;
			res[r * m + t + 2] += s02;
		}
	}
	for (; r < v; r += 3)
	{
		for (t = 0; t < h3; t++)
		{
			s00 = 0;
			s10 = 0;
			s20 = 0;
			for (q = 0; q < ah; q++)
			{
				s00 += a[r * m + q] * b[q * m + t];
				s10 += a[(r + 1) * m + q] * b[q * m + t];
				s20 += a[(r + 2) * m + q] * b[q * m + t];
			}
			res[r * m + t] += s00;
			res[(r + 1) * m + t] += s10;
			res[(r + 2) * m + t] += s20;
		}
		for (; t < h; t += 3)
		{
			s00 = 0;
			s01 = 0;
			s02 = 0;
			s10 = 0;
			s11 = 0;
			s12 = 0;
			s20 = 0;
			s21 = 0;
			s22 = 0;
			for (q = 0; q < ah; q++)
			{
				s00 += a[r * m + q] * b[q * m + t];
				s01 += a[r * m + q] * b[q * m + t + 1];
				s02 += a[r * m + q] * b[q * m + t + 2];
				s10 += a[(r + 1) * m + q] * b[q * m + t];
				s11 += a[(r + 1) * m + q] * b[q * m + t + 1];
				s12 += a[(r + 1) * m + q] * b[q * m + t + 2];
				s20 += a[(r + 2) * m + q] * b[q * m + t];
				s21 += a[(r + 2) * m + q] * b[q * m + t + 1];
				s22 += a[(r + 2) * m + q] * b[q * m + t + 2];
			}
			res[r * m + t] += s00;
			res[r * m + t + 1] += s01;
			res[r * m + t + 2] += s02;
			res[(r + 1) * m + t] += s10;
			res[(r + 1) * m + t + 1] += s11;
			res[(r + 1) * m + t + 2] += s12;
			res[(r + 2) * m + t] += s20;
			res[(r + 2) * m + t + 1] += s21;
			res[(r + 2) * m + t + 2] += s22;
		}
	}
}

void matrix_minus(double *a, double *b, int i_block, int j_block, int k, int m, int l)
{
	int h = (i_block < k ? m : l);
	int v = (j_block < k ? m : l);
	int tmp;
	for (int i = 0; i < h; i++)
	{
		tmp = i * m;
		for (int j = 0; j < v; j++)
			a[tmp + j] = a[tmp + j] - b[tmp + j];
	}
}
void matrix_subtraction(double *matrix_1, double *unit_matrix, int n)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			matrix_1[i * n + j] -= unit_matrix[i * n + j];
}

void residual(double &r1, double &r2, double *matrix, double *inverse_matrix,
							double *result_matrix, double *unit_matrix, int n, int m, int result_calc, int r)
{
	if (result_calc == SUCCESS) //  Если алгоритм применим
	{
		if (inverse_matrix)
		{
			printf("inversed matrix\n");
			print_matrix(r, n, inverse_matrix);
			if (n <= 11000)
			{
				mult_matrix(matrix, inverse_matrix, result_matrix, n, m); // A * A^(-1)
				matrix_subtraction(result_matrix, unit_matrix, n);
				r1 = norma(result_matrix, n); // ||A * A^(-1) - E ||

				mult_matrix(inverse_matrix, matrix, result_matrix, n, m); // A^(-1) * A
				matrix_subtraction(result_matrix, unit_matrix, n);
				r2 = norma(result_matrix, n); // ||A * A^(-1) - E ||
			}
			else
			{
				r1 = 0;
				r2 = 0;
			}
		}
	}
	else
	{
		r1 = -1;
		r2 = -1;
	}
}

double
find_residual(double *a, double *b, double *block, double *sum_array, int n, int m)
{
	int i = 0, p = 0, q = 0;
	int v = 0, h = 0;
	int r = 0, t = 0;
	int s = 0;
	int ah = 0;
	int z = 0;
	int v3 = 0, h3 = 0;
	double res = 0;
	double s00 = 0, s01 = 0, s02 = 0, s10 = 0, s20 = 0, s11 = 0, s12 = 0,
				 s21 = 0, s22 = 0, sum = 0;
	double *pa = nullptr, *pb = nullptr;
	int k = n / m;
	int l = n % m;
	int bl = (l != 0 ? k + 1 : k);

	for (int j = 0; j < bl; j++)
	{
		h = (j < k ? m : l);
		for (p = 0; p < h; p++)
			sum_array[p] = 0;

		for (i = 0; i < bl; i++)
		{
			v = (i < k ? m : l);
			for (p = 0; p < m * m; p++)
				block[p] = 0;

			for (s = 0; s < bl; s++)
			{
				ah = (s < k) ? m : l;
				pa = a + i * n * m + s * m;
				pb = b + s * n * m + j * m;
				v3 = v % 3;
				h3 = h % 3;
				for (r = 0; r < v3; r++)
				// первые не поделившиеся на 3 строки обслуживаем
				{
					for (t = 0; t < h3; t++)
					{
						sum = 0;
						for (q = 0; q < ah; q++)
							sum += pa[r * n + q] * pb[q * n + t];
						block[r * m + t] += sum;
					}
					for (; t < h; t += 3)
					// Внутри строки сделали кусок, который не
					// поделился нацело |****|_______________|
					{
						s00 = 0;
						s01 = 0;
						s02 = 0;
						for (q = 0; q < ah; q++)
						{
							s00 += pa[r * n + q] * pb[q * n + t];
							s01 += pa[r * n + q] * pb[q * n + t + 1];
							s02 += pa[r * n + q] * pb[q * n + t + 2];
						}
						block[r * m + t] += s00;
						block[r * m + t + 1] += s01;
						block[r * m + t + 2] += s02;
					}
				}
				for (; r < v; r += 3)
				{
					for (t = 0; t < h3; t++)
					{
						s00 = 0;
						s10 = 0;
						s20 = 0;
						for (q = 0; q < ah; q++)
						{
							s00 += pa[r * n + q] * pb[q * n + t];
							s10 += pa[(r + 1) * n + q] * pb[q * n + t];
							s20 += pa[(r + 2) * n + q] * pb[q * n + t];
						}
						block[r * m + t] += s00;
						block[(r + 1) * m + t] += s10;
						block[(r + 2) * m + t] += s20;
					}
					for (; t < h; t += 3) // здесь все делится на 3
					{
						s00 = 0;
						s01 = 0;
						s02 = 0;
						s10 = 0;
						s11 = 0;
						s12 = 0;
						s20 = 0;
						s21 = 0;
						s22 = 0; // 9 регистров
						for (q = 0; q < ah; q++)
						{
							s00 += pa[r * n + q] * pb[q * n + t];
							s01 += pa[r * n + q] * pb[q * n + t + 1];
							s02 += pa[r * n + q] * pb[q * n + t + 2];
							s10 += pa[(r + 1) * n + q] * pb[q * n + t];
							s11 += pa[(r + 1) * n + q] * pb[q * n + t + 1];
							s12 += pa[(r + 1) * n + q] * pb[q * n + t + 2];
							s20 += pa[(r + 2) * n + q] * pb[q * n + t];
							s21 += pa[(r + 2) * n + q] * pb[q * n + t + 1];
							s22 += pa[(r + 2) * n + q] * pb[q * n + t + 2];
						}
						block[r * m + t] += s00;
						block[r * m + t + 1] += s01;
						block[r * m + t + 2] += s02;
						block[(r + 1) * m + t] += s10;
						block[(r + 1) * m + t + 1] += s11;
						block[(r + 1) * m + t + 2] += s12;
						block[(r + 2) * m + t] += s20;
						block[(r + 2) * m + t + 1] += s21;
						block[(r + 2) * m + t + 2] += s22;
					}
				}
			}
			if (i == j)
				for (p = 0; p < v; p++)
					block[p * m + p] -= 1;

			for (p = 0; p < h; p++)
				for (z = 0; z < v; z++)
					sum_array[p] += fabs(block[z * m + p]);
		}

		for (i = 0; i < h; i++)
			if (sum_array[i] > res)
				res = sum_array[i];
	}
	return res;
}

int gauss_method(double *matrix, double *inversed_matrix, int *index, int *index_block,
								 double *block1, double *invert_block, double *block3, double *block4, int n,
								 int m, double matrix_norm)
{
	int i, j;
	int k = n / m;
	int l = n % m;
	int index_min_block = 0, count_blocks = 0, tmp = 0;
	int index_0 = 0, index_1 = 0;
	double norm_min_block = 0, tmp_norm_block = 0;
	count_blocks = (l > 0 ? k + 1 : k);

	for (i = 0; i < k; i++)
	{
		index_min_block = -1;
		norm_min_block = 0;

		for (j = i; j < k; j++) // Нахождение блока с наименьшней нормой
		{
			tmp_norm_block = 0;
			get_block(matrix, block1, i, j, n, m);
			if (gauss_classic_row(block1, invert_block, index, m, matrix_norm, m) != IRREVERSIBLE)
			{
				tmp_norm_block = norma(invert_block, m);
				if (norm_min_block > tmp_norm_block || index_min_block == -1)
				{
					norm_min_block = tmp_norm_block;
					index_min_block = j;
				}
			}
		}

		if (index_min_block == -1) // Если все блоки вырожденные, метод неприменим
			return CANNOT_SOLVE;

		for (j = 0; j < count_blocks; j++) // Переставляем блоки
		{
			get_block(matrix, block1, j, i, n, m);
			get_block(matrix, invert_block, j, index_min_block, n, m);
			set_block(matrix, block1, j, index_min_block, n, m);
			set_block(matrix, invert_block, j, i, n, m);
		}

		for (j = 0; j < m; j++)
		{
			tmp = index_block[index_min_block * m + j];
			index_block[index_min_block * m + j] = index_block[i * m + j];
			index_block[i * m + j] = tmp;
		}
		get_block(matrix, block1, i, i, n, m);
		gauss_classic_row(block1, invert_block, index, m, matrix_norm, m);
		get_block(matrix, block1, i, i, n, m);
		set_unit_zero_block(matrix, block1, i, i, n, m);

		for (int j = i + 1; j < count_blocks; j++) // Умножаем i строку на обратный элемент
		{
			get_block(matrix, block1, i, j, n, m);
			j == k ? mult(invert_block, block1, block3, m, m, l, m) : mult(invert_block, block1, block3, m, m, m, m);
			set_block(matrix, block3, i, j, n, m);
		}

		for (int j = 0; j < count_blocks; j++)
		{
			get_block(inversed_matrix, block1, i, j, n, m); // В присоединенной матрице
			j == k ? mult(invert_block, block1, block3, m, m, l, m) : mult(invert_block, block1, block3, m, m, m, m);
			set_block(inversed_matrix, block3, i, j, n, m);
		}

		for (int r = i + 1; r < count_blocks; r++)
		{
			get_block(matrix, block1, r, i, n, m);
			get_block(matrix, block3, r, i, n, m);
			set_unit_zero_block(matrix, block3, r, i, n, m);

			for (int j = 0; j < count_blocks; j++)
			{
				get_block(matrix, invert_block, i, j, n, m);
				if (r == k)
					j == k ? mult(block1, invert_block, block4, l, m, l, m) : mult(block1, invert_block, block4, l, m, m, m);
				else
					j == k ? mult(block1, invert_block, block4, m, m, l, m) : mult(block1, invert_block, block4, m, m, m, m);

				get_block(matrix, block3, r, j, n, m);

				matrix_minus(block3, block4, r, j, k, m, l);
				set_block(matrix, block3, r, j, n, m);
			}

			for (int j = 0; j < count_blocks; j++)
			{
				get_block(inversed_matrix, invert_block, i, j, n, m);
				if (r == k)
					j == k ? mult(block1, invert_block, block4, l, m, l, m) : mult(block1, invert_block, block4, l, m, m, m);
				else
					j == k ? mult(block1, invert_block, block4, m, m, l, m) : mult(block1, invert_block, block4, m, m, m, m);

				get_block(inversed_matrix, block3, r, j, n, m);
				matrix_minus(block3, block4, r, j, k, m, l);
				set_block(inversed_matrix, block3, r, j, n, m);
			}
		}
	}

	if (l != 0)
	{
		get_block(matrix, block1, k, k, n, m);
		if (gauss_classic_row(block1, invert_block, index, l, matrix_norm, m) == IRREVERSIBLE)
			return CANNOT_SOLVE;
		set_unit_zero_block(matrix, block4, k, k, n, m);
		for (int j = 0; j <= k; j++)
		{
			if (j == k)
			{
				get_block(inversed_matrix, block1, k, j, n, m);
				mult(invert_block, block1, block3, l, l, l, m);
				set_block(inversed_matrix, block3, k, j, n, m);
			}
			else
			{
				get_block(inversed_matrix, block1, k, j, n, m);
				mult(invert_block, block1, block3, l, l, m, m);
				set_block(inversed_matrix, block3, k, j, n, m);
			}
		}
	}

	for (i = n - 1; i >= 0; i--)
		for (k = i - 1; k >= 0; k--)
		{
			index_0 = k * n;
			index_1 = i * n;
			for (j = 0; j < n; j++)
			{
				inversed_matrix[index_0 + j] -= inversed_matrix[index_1 + j] * matrix[index_0 + i];
			}
			matrix[index_0 + i] = 0;
		}

	for (int i = 0; i < n; i++)
	{
		index_0 = index_block[i] * n;
		index_1 = i * n;
		for (int j = 0; j < n; j++)
			matrix[index_0 + j] = inversed_matrix[index_1 + j];
	}

	for (int i = 0; i < n; i++)
	{
		index_0 = i * n;
		for (int j = 0; j < n; j++)
			inversed_matrix[index_0 + j] = matrix[index_0 + j];
	}
	return SUCCESS;
}

void set_unit_zero_block(double *matrix, double *block, int i_block, int j_block, int n, int m)
{
	int tmp;
	// делаем единичным диагональный блок и нулевой под диагональным
	for (int i = 0; i < m; i++)
	{
		tmp = i * m;
		for (int j = 0; j < m; j++)
			j == i &&i_block == j_block ? block[tmp + j] = 1 : block[tmp + j] = 0;
	}
	set_block(matrix, block, i_block, j_block, n, m);
}

int gauss_classic_row(double *matrix, double *inverse_matrix, int *index, int n, double matrix_norm, int row_ind)
{
	int max_index = 0;
	double tmp_ = 0, max = 0;
	int index_0 = 0;
	int index_1 = 0;
	int swap = 0;

	for (int i = 0; i < n; i++) // Присоединенная матрица
	{
		index_0 = i * row_ind;
		for (int j = 0; j < n; j++)
			i == j ? inverse_matrix[index_0 + j] = 1. : inverse_matrix[index_0 + j] = 0.;
	}

	for (int i = 0; i < n; i++)
		index[i] = i;

	for (int i = 0; i < n; i++) // Прямой ход метода Гаусса
	{
		index_0 = i * row_ind;
		max = fabs(matrix[index_0 + i]); // нахождение главного элемента
		max_index = i;

		for (int j = i + 1; j < n; j++)
			if (max < fabs(matrix[index_0 + j]))
			{
				max = fabs(matrix[index_0 + j]);
				max_index = j;
			}

		swap = index[i];
		index[i] = index[max_index];
		index[max_index] = swap;

		for (int j = 0; j < n; j++) // переставляем столбцы местами
		{
			index_1 = j * row_ind;
			tmp_ = matrix[index_1 + i];
			matrix[index_1 + i] = matrix[index_1 + max_index];
			matrix[index_1 + max_index] = tmp_;
		}

		if (fabs(matrix[index_0 + i]) < matrix_norm * EPS) // если элемент нулевой, метод неприменим
		{
			return IRREVERSIBLE;
		}

		// tmp_ = matrix[index_0 + i];
		tmp_ = 1 / matrix[index_0 + i];
		matrix[index_0 + i] = 1.0;
		for (int j = i + 1; j < n; j++) // умножаем i строку на обратный элемент
			matrix[index_0 + j] *= tmp_;

		for (int j = 0; j < n; j++) // присоединенная матрица. умножаем i строку на обратный элемент
			inverse_matrix[index_0 + j] *= tmp_;

		for (int j = i + 1; j < n; j++)
		{
			tmp_ = matrix[j * row_ind + i];
			for (int k = i; k < n; k++) // вычитание строки
				matrix[j * row_ind + k] -= matrix[index_0 + k] * tmp_;
			for (int k = 0; k < n; k++)
				inverse_matrix[j * row_ind + k] -= inverse_matrix[index_0 + k] * tmp_; // присоединенная матрица. вычитание строки умноженной на число
		}
	}

	for (int i = 0; i < n; i++) // Обратный ход
		for (int j = n - 1; j >= 0; j--)
		{
			index_0 = j * row_ind;
			tmp_ = inverse_matrix[index_0 + i];
			for (int k = j + 1; k < n; k++)
				tmp_ -= matrix[index_0 + k] * inverse_matrix[k * row_ind + i];
			inverse_matrix[index_0 + i] = tmp_;
		}

	for (int i = 0; i < n; i++)
	{
		index_0 = i * row_ind;
		index_1 = index[i] * row_ind;
		for (int j = 0; j < n; j++)
			matrix[index_1 + j] = inverse_matrix[index_0 + j];
	}

	for (int i = 0; i < n; i++)
	{
		index_0 = i * row_ind;
		for (int j = 0; j < n; j++)
			inverse_matrix[index_0 + j] = matrix[index_0 + j];
	}

	return SUCCESS;
}

void mult_matrix(double *a, double *b, double *c, int n, int m)
{
	int k = 0, l = 0, i = 0, j = 0;
	int v = 0, h = 0;
	int r = 0, t = 0;
	int s = 0;
	int ah = 0;
	int q = 0;
	int v3 = 0, h3 = 0;
	double s00 = 0, s01 = 0, s02 = 0, s10 = 0, s20 = 0, s11 = 0, s12 = 0,
				 s21 = 0, s22 = 0, sum = 0;
	double *pa = nullptr, *pb = nullptr, *pc = nullptr;
	k = n / m;
	l = n % m;
	int bl = (l != 0 ? k + 1 : k);
	for (i = 0; i < bl; i++)
	{
		for (j = 0; j < bl; j++) // c_ij = \sum s=1^l a_is*b_sj;
		{
			// Размер c_ij v x h высота на ширину v строк , h столбцов
			v = (i < k ? m : l);
			h = (j < k ? m : l);
			pc = c + i * n * m + j * m; //верхний левый угол блока
			for (r = 0; r < v; r++)			// Весь блок матрицы С попадает в кэш
																	// память.Обнулим матрицу C.
				for (t = 0; t < h; t++)
					pc[r * n + t] = 0; // может быть pc[r * h + t] = 0;
			// Начнем цикл по S.
			for (s = 0; s < bl; s++)
			{
				ah = (s < k ? m : l);				// столбцов в A
				pa = a + i * n * m + s * m; // A_is | v x ah
				pb = b + s * n * m + j * m; // B_sj | ah x h
				// Умножаем блок A_is на B_sj Обслуживаем случаи:
				v3 = v % 3;
				h3 = h % 3;
				for (r = 0; r < v3;
						 r++) // первые не поделившиеся на 3 строки обслуживаем
				{
					for (t = 0; t < h3;
							 t++) // Внутри строки сделали кусок, который не
										// поделился нацело |****|_______________|
					{
						sum = 0; // s=0;
						for (q = 0; q < ah; q++)
							sum += pa[r * n + q] * pb[q * n + t]; // s
						pc[r * n + t] += sum;										// s
					}
					for (; t < h;
							 t += 3) // тут поделилось. Делаем разворачивание в одной
											 // строке |____|************|
					{
						s00 = 0;
						s01 = 0;
						s02 = 0;
						for (q = 0; q < ah; q++)
						{
							s00 += pa[r * n + q] * pb[q * n + t];
							s01 += pa[r * n + q] * pb[q * n + (t + 1)];
							s02 += pa[r * n + q] * pb[q * n + (t + 2)];
						}
						pc[r * n + t] += s00;
						pc[r * n + (t + 1)] += s01;
						pc[r * n + (t + 2)] += s02;
					}
				}
				for (; r < v; r += 3) // обслуживаем другие случаи.
				{
					for (t = 0; t < h3; t++)
					{
						s00 = 0;
						s10 = 0;
						s20 = 0;
						for (q = 0; q < ah; q++)
						{
							s00 += pa[r * n + q] * pb[q * n + t];
							s10 += pa[(r + 1) * n + q] * pb[q * n + t];
							s20 += pa[(r + 2) * n + q] * pb[q * n + t];
						}
						pc[r * n + t] += s00;
						pc[(r + 1) * n + t] += s10;
						pc[(r + 2) * n + t] += s20;
					}
					for (; t < h; t += 3) // здесь все делится на 3
					{
						s00 = 0;
						s01 = 0;
						s02 = 0;
						s10 = 0;
						s11 = 0;
						s12 = 0;
						s20 = 0;
						s21 = 0;
						s22 = 0; // 9 регистров
						for (q = 0; q < ah; q++)
						{
							s00 += pa[r * n + q] * pb[q * n + t];
							s01 += pa[r * n + q] * pb[q * n + (t + 1)];
							s02 += pa[r * n + q] * pb[q * n + (t + 2)];
							s10 += pa[(r + 1) * n + q] * pb[q * n + t];
							s11 += pa[(r + 1) * n + q] * pb[q * n + (t + 1)];
							s12 += pa[(r + 1) * n + q] * pb[q * n + (t + 2)];
							s20 += pa[(r + 2) * n + q] * pb[q * n + t];
							s21 += pa[(r + 2) * n + q] * pb[q * n + (t + 1)];
							s22 += pa[(r + 2) * n + q] * pb[q * n + (t + 2)];
						}
						pc[r * n + t] += s00;
						pc[r * n + (t + 1)] += s01;
						pc[r * n + (t + 2)] += s02;
						pc[(r + 1) * n + t] += s10;
						pc[(r + 1) * n + (t + 1)] += s11;
						pc[(r + 1) * n + (t + 2)] += s12;
						pc[(r + 2) * n + t] += s20;
						pc[(r + 2) * n + (t + 1)] += s21;
						pc[(r + 2) * n + (t + 2)] += s22;
					}
				}
			}
		}
	}
}