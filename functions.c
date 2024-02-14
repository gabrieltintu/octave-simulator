// Tintu Gabriel-Claudiu Grupa 323CA 2023 - 2024

#include "functions.h"

// The remainder when divided by 10007
int mod(int n)
{
	n = ((n % 10007) + 10007) % 10007;
	return n;
}

// Allocating vectors for the number of rows and columns of matrices
void allocate_rows_collumns(int **vek_n, int **vek_m, int cnt)
{
	(*vek_n) = realloc((*vek_n), cnt * sizeof(int));
	if (!(*vek_n)) {
		printf("Could not allocate.\n");
		exit(-1);
	}
	(*vek_m) = realloc((*vek_m), cnt * sizeof(int));
	if (!(*vek_m)) {
		printf("Could not allocate.\n");
		exit(-1);
	}
}

// Allocate memory and read matrices
void allocate_list(int ****a, int n, int m, int cnt)
{
	(*a) = realloc((*a), cnt * sizeof(int **));
	if (!(*a)) {
		printf("Could not allocate.\n");
		exit(-1);
	}
	(*a)[cnt - 1] = malloc(n * sizeof(int *));
	if (!(*a)[cnt - 1]) {
		printf("Could not allocate.\n");
		exit(-1);
	}
	for (int i = 0; i < n; i++) {
		(*a)[cnt - 1][i] = malloc(m * sizeof(int));
		if (!(*a)[cnt - 1][i]) {
			printf("Could not allocate.\n");
			exit(-1);
		}
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++)
			scanf("%d", &(*a)[cnt - 1][i][j]);
	}
}

// Allocate memory and read the rows and columns vectors for redimension
void allocate_rows_colls_red(int **vek_l, int **vek_c, int *nlin, int *ncol)
{
	scanf("%d", &(*nlin));
	// vek_l rows indices
	(*vek_l) = malloc((*nlin) * sizeof(int));
	if (!(*vek_l)) {
		printf("Could not allocate.\n");
		exit(-1);
	}
	for (int i = 0; i < (*nlin); i++)
		scanf("%d", &(*vek_l)[i]);
	scanf("%d", &(*ncol));
	// vek_c collumns indices
	(*vek_c) = malloc((*ncol) * sizeof(int));
	if (!(*vek_c)) {
		printf("Could not allocate.\n");
		exit(-1);
	}
	for (int i = 0; i < (*ncol); i++)
		scanf("%d", &(*vek_c)[i]);
}

// Print matrices
void print_matrices(int x, int ***a, int *vek_n, int *vek_m, int cnt)
{
	if (x < cnt && x >= 0) {
		for (int i = 0; i < vek_n[x]; i++) {
			for (int j = 0; j < vek_m[x]; j++)
				printf("%d ", a[x][i][j]);
			printf("\n");
		}
	} else {
		printf("No matrix with the given index\n");
	}
}

// Allocate an auxiliary matrix with n rows and m columns to store values
// that will be copied into the matrix with index x
// Using realloc to resize matrix x and then copying
// values from the auxiliary matrix, followed by freeing it
void red_mat(int ****a, int x, int n, int m, int *vek_l, int *vek_c, int nr)
{
	int **mat_aux = malloc(n * sizeof(int *));
	if (!mat_aux) {
		printf("Could not allocate.\n");
		exit(-1);
	}

	for (int i = 0; i < n; i++) {
		mat_aux[i] = malloc(m * sizeof(int));
		if (!mat_aux[i]) {
			printf("Could not allocate.\n");
			exit(-1);
		}
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++)
			mat_aux[i][j] = (*a)[x][vek_l[i]][vek_c[j]];
	}
	for (int i = 0; i < nr; i++)
		free((*a)[x][i]);
	free((*a)[x]);
	(*a)[x] = mat_aux;
}

// The function multiplies the two matrices with indices x and y
// into an auxiliary matrix, then allocates additional memory
// for the product matrix at the end of the matrix list
// The values from the auxiliary matrix are copied into the product matrix
void multiply_mat(int ****a, int x, int y, int **vek_n, int **vek_m, int cnt)
{
	// n1 = nr of rows for the matrix qith index x
	// n2 = nr of rows for the matrix with index y
	// same for m1, m2 (nr of columns)
	int n1 = (*vek_n)[x], m1 = (*vek_m)[x], m2 = (*vek_m)[y];
	int n = n1;
	int m = m2;
	int **mat_aux = malloc(n * sizeof(int *));
	if (!mat_aux) {
		printf("Could not allocate.\n");
		exit(-1);
	}

	for (int i = 0; i < n; i++) {
		mat_aux[i] = malloc(m * sizeof(int));
		if (!mat_aux[i]) {
			printf("Could not allocate.\n");
			exit(-1);
		}
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++)
			mat_aux[i][j] = 0;
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			for (int k = 0; k < m1; k++) {
				mat_aux[i][j] += (*a)[x][i][k] * (*a)[y][k][j];
				mat_aux[i][j] = mod(mat_aux[i][j]);
			}
		}
	}
	(*a) = realloc((*a), cnt * sizeof(int **));
	if (!(*a)) {
		printf("Could not allocate.\n");
		exit(-1);
	}
	(*a)[cnt - 1] = malloc(n * sizeof(int *));
	if (!(*a)[cnt - 1]) {
		printf("Could not allocate.\n");
		exit(-1);
	}
	for (int i = 0; i < n; i++) {
		(*a)[cnt - 1][i] = malloc(m * sizeof(int));
		if (!(*a)[cnt - 1][i]) {
			printf("Could not allocate.\n");
			exit(-1);
		}
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++)
			(*a)[cnt - 1][i][j] = mat_aux[i][j];
	}

	// Allocate additional space at the end of the vectors
	// that store the number of columns and rows
	// Adding the number of rows and columns of the product matrix
	// to the corresponding vectors

	allocate_rows_collumns(&(*vek_n), &(*vek_m), cnt);
	(*vek_n)[cnt - 1] = n;
	(*vek_m)[cnt - 1] = m;

	for (int i = 0; i < n; i++)
		free(mat_aux[i]);
	free(mat_aux);
}

// Calculate the sums of matrices in the list
// and storing them in the v_sum vector
void matrix_sum(int ***a, int **v_sum, int *vek_n, int *vek_m, int cnt)
{
	(*v_sum) = realloc((*v_sum), cnt * sizeof(int));
	if (!(*v_sum)) {
		printf("Could not allocate.\n");
		exit(-1);
	}
	for (int i = 0; i < cnt; i++)
		(*v_sum)[i] = 0;
	for (int k = 0; k < cnt; k++)  {
		for (int i = 0; i < vek_n[k]; i++) {
			for (int j = 0; j < vek_m[k]; j++)
				(*v_sum)[k] += mod(a[k][i][j]);
		}
	}
	for (int i = 0; i < cnt; i++)
		(*v_sum)[i] = mod((*v_sum)[i]);
}

// Sort matrices, vectors of sums, rows and columns
void sort_matrices(int ****a, int **v_sum, int **vek_n, int **vek_m, int cnt)
{
	for (int k = 0; k < cnt - 1; k++) {
		for (int l = k + 1; l < cnt; l++) {
			if ((*v_sum)[k] > (*v_sum)[l]) {
				int aux;
				aux = (*v_sum)[k];
				(*v_sum)[k] = (*v_sum)[l];
				(*v_sum)[l] = aux;
				int **mat_aux = (*a)[k];
				(*a)[k] = (*a)[l];
				(*a)[l] = mat_aux;
				int naux = (*vek_n)[k];
				(*vek_n)[k] = (*vek_n)[l];
				(*vek_n)[l] = naux;
				int maux = (*vek_m)[k];
				(*vek_m)[k] = (*vek_m)[l];
				(*vek_m)[l] = maux;
			}
		}
	}
}

// Store the transposed matrix with index x in an auxiliary matrix
// Reallocate memory space to copy the transpose into the list
void transpose_matrix(int ****a, int **vek_n, int **vek_m, int x)
{
	int n = (*vek_n)[x];
	int m = (*vek_m)[x];
	int **mat_aux = malloc(m * sizeof(int *));
	if (!mat_aux) {
		printf("Could not allocate.\n");
		exit(-1);
	}

	for (int i = 0; i < m; i++) {
		mat_aux[i] = malloc(n * sizeof(int));
		if (!mat_aux[i]) {
			printf("Could not allocate.\n");
			exit(-1);
		}
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++)
			mat_aux[j][i] = (*a)[x][i][j];
	}
	for (int i = 0; i < n; i++)
		free((*a)[x][i]);
	free((*a)[x]);

	(*a)[x] = mat_aux;

	(*vek_n)[x] = m;
	(*vek_m)[x] = n;
}

// Calculate the product of two matrices and return it
int **mat_product(int **mat_aux1, int **mat_aux2, int n, int ***rez)
{
	int **mat_aux = malloc(n * sizeof(int *));
	if (!mat_aux) {
		printf("Could not allocate.\n");
		exit(-1);
	}

	for (int i = 0; i < n; i++) {
		mat_aux[i] = calloc(n, sizeof(int));
		if (!mat_aux[i]) {
			printf("Could not allocate.\n");
			exit(-1);
		}
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				mat_aux[i][j] += mat_aux1[i][k] * mat_aux2[k][j];
				mat_aux[i][j] = mod(mat_aux[i][j]);
			}
		}
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			(*rez)[i][j] = mat_aux[i][j];
	}
	for (int i = 0; i < n; i++)
		free(mat_aux[i]);
	free(mat_aux);
	return *rez;
}

// Raise the matrix with the index x to the power p
int **power_matrix(int ***a, int n, int p, int x, int ***rez)
{
	int **z = a[x];
	if (p == 1)
		return z;
	if (p % 2 == 1) {
		int **mat_aux = power_matrix(a, n, (p - 1) / 2, x, rez);
		return mat_product(z, mat_product(mat_aux, mat_aux, n, rez), n, rez);
	}
	int **mat_aux = power_matrix(a, n, p / 2, x, rez);
	return mat_product(mat_aux, mat_aux, n, rez);
}

// Delete from the list and free the matrix with the index x
void free_matrix(int ****a, int **vek_n, int **vek_m, int x, int *cnt)
{
	for (int i = x; i < (*cnt) - 1; i++) {
		if (i == x) {
			for (int j = 0; j < (*vek_n)[i]; j++)
				free((*a)[i][j]);
			free((*a)[i]);
		}
		(*a)[i] = (*a)[i + 1];
	}
	for (int i = x; i < (*cnt) - 1; i++) {
		(*vek_n)[i] = (*vek_n)[i + 1];
		(*vek_m)[i] = (*vek_m)[i + 1];
	}
	(*cnt)--;
	if (x == (*cnt)) {
		for (int i = 0; i < (*vek_n)[(*cnt)]; i++)
			free((*a)[(*cnt)][i]);
		free((*a)[(*cnt)]);
	}
	(*a) = realloc((*a), (*cnt) * sizeof(int **));
	(*vek_n) = realloc((*vek_n), (*cnt) * sizeof(int));
	(*vek_m) = realloc((*vek_m), (*cnt) * sizeof(int));
}

// Deallocate the matrix list, vectors for rows and columns
void deallocate(int ****a, int **vek_n, int **vek_m, int cnt, int **v_sum)
{
	for (int i = 0; i < cnt; i++) {
		for (int j = 0; j < (*vek_n)[i]; j++)
			free((*a)[i][j]);
		free((*a)[i]);
	}
	free((*a));
	free((*vek_n));
	free((*vek_m));
	free((*v_sum));
}

// Conditions to display messages for the R operation
void verify_R(int x, int p, int *vek_n, int *vek_m, int cnt, int ****a)
{
	if (x >= cnt || x < 0) {
		printf("No matrix with the given index\n");
	} else if (p < 0) {
		printf("Power should be positive\n");
	} else if (vek_n[x] != vek_m[x]) {
		printf("Cannot perform matrix multiplication\n");
	} else if (x < cnt && x >= 0 && p > 0 && vek_n[x] == vek_m[x]) {
		int nr = vek_n[x];
		int **rez = malloc(nr * sizeof(int *));
		int **aux = (*a)[x];
		for (int i = 0; i < nr; i++)
			rez[i] = calloc(nr, sizeof(int));
		power_matrix((*a), nr, p, x, &rez);
		(*a)[x] = rez;
		for (int i = 0; i < nr; i++)
			free(aux[i]);
		free(aux);
	}
}

// Allocate memory for a n x n matrix
int **allocate_matrix(int n)
{
	int **aa = malloc(n * sizeof(int *));
	if (!aa) {
		printf("Could not allocate.\n");
		exit(-1);
	}
	for (int i = 0; i < n ; i++) {
		aa[i] = malloc(n * sizeof(int));
		if (!aa[i]) {
			printf("Could not allocate.\n");
			exit(-1);
		}
	}
	return aa;
}

// Add 2 matrices
int **add(int **aux1, int **aux2, int n)
{
	int **aux3 = allocate_matrix(n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			aux3[i][j] = aux1[i][j] + aux2[i][j];
	return aux3;
}

// Substract 2 matrices
int **substract(int **aux1, int **aux2, int n)
{
	int **aux3 = allocate_matrix(n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			aux3[i][j] = aux1[i][j] - aux2[i][j];
	return aux3;
}

// Free matrices used for the strassen2
void free_strassen(int **m1, int **m2, int **m3, int **m4, int n)
{
	for (int i = 0; i < n; i++)
		free(m1[i]);
	free(m1);
	for (int i = 0; i < n; i++)
		free(m2[i]);
	free(m2);
	for (int i = 0; i < n; i++)
		free(m3[i]);
	free(m3);
	for (int i = 0; i < n; i++)
		free(m4[i]);
	free(m4);
}

// Functions values_1,2,3,4 will copy the values from the matrices
// to multiply, in those 8 smaller matrices
void values_1(int ***sa1, int ***sb1, int **aux1, int **aux2, int n)
{
	(*sa1) = allocate_matrix(n / 2);
	(*sb1) = allocate_matrix(n / 2);
	for (int i = 0; i < n / 2; i++) {
		for (int j = 0; j < n / 2; j++) {
			(*sa1)[i][j] = aux1[i][j];
			(*sb1)[i][j] = aux2[i][j];
		}
	}
}

void values_2(int ***sa2, int ***sb2, int **aux1, int **aux2, int n)
{
	(*sa2) = allocate_matrix(n / 2);
	(*sb2) = allocate_matrix(n / 2);
	for (int i = 0; i < n / 2; i++) {
		for (int j = n / 2; j < n; j++) {
			(*sa2)[i][j - n / 2] = aux1[i][j];
			(*sb2)[i][j - n / 2] = aux2[i][j];
		}
	}
}

void values_3(int ***sa3, int ***sb3, int **aux1, int **aux2, int n)
{
	(*sa3) = allocate_matrix(n / 2);
	(*sb3) = allocate_matrix(n / 2);
	for (int i = n / 2; i < n; i++) {
		for (int j = 0; j < n / 2; j++) {
			(*sa3)[i - n / 2][j] = aux1[i][j];
			(*sb3)[i - n / 2][j] = aux2[i][j];
		}
	}
}

void values_4(int ***sa4, int ***sb4, int **aux1, int **aux2, int n)
{
	(*sa4) = allocate_matrix(n / 2);
	(*sb4) = allocate_matrix(n / 2);
	for (int i = n / 2; i < n; i++) {
		for (int j = n / 2; j < n; j++) {
			(*sa4)[i - n / 2][j - n / 2] = aux1[i][j];
			(*sb4)[i - n / 2][j - n / 2] = aux2[i][j];
		}
	}
}

// Place the values calculated using the Strassen algorithm into the matrix mat
int **val_mat(int **mat, int **cm1, int **cm2, int **cm3, int **cm4, int n)
{
	for (int i = 0; i < n / 2; i++) {
		for (int j = 0; j < n / 2; j++)
			mat[i][j] = mod(cm1[i][j]);
	}
	for (int i = 0; i < n / 2; i++) {
		for (int j = n / 2; j < n; j++)
			mat[i][j] = mod(cm2[i][j - n / 2]);
	}
	for (int i = n / 2; i < n; i++) {
		for (int j = 0; j < n / 2; j++)
			mat[i][j] = mod(cm3[i - n / 2][j]);
	}
	for (int i = n / 2; i < n; i++) {
		for (int j = n / 2; j < n; j++)
			mat[i][j] = mod(cm4[i - n / 2][j - n / 2]);
	}
	return mat;
}

// Multiply two matrices using the Strassen algorithm
// If n = 2, the algorithm will effectively calculate using the formula
// Otherwise, if n > 2, the algorithm will split the two matrices
// into four matrices of size n / 2, and operations (additions, subtractions)
// will be applied
// For multiplications, the strassen2 function will be called again
int **strassen2(int **aux1, int **aux2, int n)
{
	if (n == 2) {
		int **aux3 = allocate_matrix(n);
		int ad1 = 0, ad2 = 0;
		ad1 = aux1[0][0] + aux1[1][1];
		ad2 = aux2[0][0] + aux2[1][1];
		int m1 = mod(mod(ad1) * mod(ad2));
		int m2 = mod(mod(aux1[1][0] + aux1[1][1]) * mod(aux2[0][0]));
		int m3 = mod(mod(aux1[0][0]) * mod(aux2[0][1] - aux2[1][1]));
		int m4 = mod(mod(aux1[1][1])  * mod(aux2[1][0] - aux2[0][0]));
		int m5 = mod(mod(aux1[0][0] + aux1[0][1]) * mod(aux2[1][1]));
		ad1 = aux1[1][0] - aux1[0][0];
		ad2 = aux2[0][0] + aux2[0][1];
		int m6 = mod(mod(ad1) * mod(ad2));
		ad1 = aux1[0][1] - aux1[1][1];
		ad2 = aux2[1][0] + aux2[1][1];
		int m7 = mod(mod(ad1) * mod(ad2));
		int c1 = mod(m1 + m4 - m5 + m7);
		int c2 = mod(m3 + m5);
		int c3 = mod(m2 + m4);
		int c4 = mod(m1 - m2 + m3 + m6);
		aux3[0][0] = mod(c1);
		aux3[0][1] = mod(c2);
		aux3[1][0] = mod(c3);
		aux3[1][1] = mod(c4);
		return aux3;
	}
	int **sa1, **sa2, **sa3, **sa4, **sb1, **sb2, **sb3, **sb4;

	values_1(&sa1, &sb1, aux1, aux2, n);
	values_2(&sa2, &sb2, aux1, aux2, n);
	values_3(&sa3, &sb3, aux1, aux2, n);
	values_4(&sa4, &sb4, aux1, aux2, n);

	int **maux1 = add(sa1, sa4, n / 2);
	int **maux2 = add(sb1, sb4, n / 2);
	int **maux3 = add(sa3, sa4, n / 2);
	int **maux4 = substract(sb2, sb4, n / 2);
	int **maux5 = substract(sb3, sb1, n / 2);
	int **maux6 = add(sa1, sa2, n / 2);
	int **maux7 = substract(sa3, sa1, n / 2);
	int **maux8 = add(sb1, sb2, n / 2);
	int **maux9 = substract(sa2, sa4, n / 2);
	int **maux10 = add(sb3, sb4, n / 2);
	int **mat1 = strassen2(maux1, maux2, n / 2);
	int **mat2 = strassen2(maux3, sb1, n / 2);
	int **mat3 = strassen2(sa1, maux4, n / 2);
	int **mat4 = strassen2(sa4, maux5, n / 2);
	int **mat5 = strassen2(maux6, sb4, n / 2);
	int **mat6 = strassen2(maux7, maux8, n / 2);
	int **mat7 = strassen2(maux9, maux10, n / 2);
	int **maux11 = add(mat1, mat4, n / 2);
	int **maux12 = substract(maux11, mat5, n / 2);
	int **maux13 = substract(mat1, mat2, n / 2);
	int **maux14 = add(maux13, mat3, n / 2);
	int **cm1 = add(maux12, mat7, n / 2);
	int **cm2 = add(mat3, mat5, n / 2);
	int **cm3 = add(mat2, mat4, n / 2);
	int **cm4 = add(maux14, mat6, n / 2);
	free_strassen(maux1, maux2, maux3, maux4, n / 2);
	free_strassen(maux5, maux6, maux7, maux8, n / 2);
	free_strassen(maux9, maux10, mat1, mat2, n / 2);
	free_strassen(mat3, mat4, mat5, mat6, n / 2);
	free_strassen(mat7, maux11, maux12, maux13, n / 2);
	for (int i = 0; i < n / 2; i++)
		free(maux14[i]);
	free(maux14);

	int **mat = allocate_matrix(n);
	mat = val_mat(mat, cm1, cm2, cm3, cm4, n);

	free_strassen(cm1, cm2, cm3, cm4, n / 2);
	free_strassen(sa1, sa2, sa3, sa4, n / 2);
	free_strassen(sb1, sb2, sb3, sb4, n / 2);
	return mat;
}

// Verifies if the two matrices can be multiplied, and if yes,
// calls the strassen2 function and makes the necessary modifications
// in the matrix list and vectors for rows and columns
void verify_strassen(int ****a, int *cnt, int **vek_n, int **vek_m)
{
	int x, y;
	scanf("%d%d", &x, &y);
	int dim = (*vek_n)[x];
	if (x >= 0 && x < (*cnt) && y >= 0 && y < (*cnt)) {
		if ((*vek_m)[x] == (*vek_n)[y]) {
			int **aux1 = (*a)[x];
			int **aux2 = (*a)[y];
			int **aux3 = strassen2(aux1, aux2, dim);
			(*cnt)++;
			(*a) = realloc((*a), (*cnt) * sizeof(int **));
			allocate_rows_collumns(&(*vek_n), &(*vek_m), (*cnt));
			(*vek_n)[(*cnt) - 1] = dim;
			(*vek_m)[(*cnt) - 1] = dim;
			(*a)[(*cnt) - 1] = aux3;
		} else {
			printf("Cannot perform matrix multiplication\n");
		}
	} else {
		printf("No matrix with the given index\n");
	}
}
