// Tintu Gabriel-Claudiu 2023 - 2024

#include "functions.h"

int main(void)
{
	// cnt - number of matrices in the list
	int cnt = 0, x, y, p, n, m, nlin, ncol;
	char c;
	// a - list of matrices
	// vek_n store nr of rows for every matrix
	// vek_m store nr of columns for every matrix
	int ***a = NULL, *vek_n = NULL, *vek_m = NULL, *v_sum = NULL;
	while (1) {
		scanf("%c", &c);
		if (c == 'L') {
			cnt++;
			allocate_rows_collumns(&vek_n, &vek_m, cnt);
			scanf("%d%d", &n, &m);
			vek_n[cnt - 1] = n;
			vek_m[cnt - 1] = m;
			allocate_list(&a, n, m, cnt);
		} else if (c == 'D') {
			scanf("%d", &x);
			if (x < cnt && x >= 0)
				printf("%d %d\n", vek_n[x], vek_m[x]);
			else
				printf("No matrix with the given index\n");
		} else if (c == 'P') {
			scanf("%d", &x);
			print_matrices(x, a, vek_n, vek_m, cnt);
		} else if (c == 'C') {
			scanf("%d", &x);
			if (x < cnt && x >= 0) {
				int *vek_l, *vek_c;
				allocate_rows_colls_red(&vek_l, &vek_c, &nlin, &ncol);
				red_mat(&a, x, nlin, ncol, vek_l, vek_c, vek_n[x]);
				vek_n[x] = nlin;
				vek_m[x] = ncol;
				free(vek_l);
				free(vek_c);
			} else {
				printf("No matrix with the given index\n");
			}
		} else if (c == 'M') {
			scanf("%d%d", &x, &y);
			if (x < cnt && y < cnt && x >= 0 && y >= 0) {
				if (vek_m[x] == vek_n[y]) {
					cnt++;
					multiply_mat(&a, x, y, &vek_n, &vek_m, cnt);
				} else {
					printf("Cannot perform matrix multiplication\n");
				}
			} else {
				printf("No matrix with the given index\n");
			}
		} else if (c == 'O') {
			matrix_sum(a, &v_sum, vek_n, vek_m, cnt);
			sort_matrices(&a, &v_sum, &vek_n, &vek_m, cnt);
		} else if (c == 'T') {
			scanf("%d", &x);
			if (x < cnt && x >= 0)
				transpose_matrix(&a, &vek_n, &vek_m, x);
			else
				printf("No matrix with the given index\n");
		} else if (c == 'R') {
			scanf("%d%d", &x, &p);
			verify_R(x, p, vek_n, vek_m, cnt, &a);
		} else if (c == 'F') {
			scanf("%d", &x);
			if (x < cnt && x >= 0)
				free_matrix(&a, &vek_n, &vek_m, x, &cnt);
			else
				printf("No matrix with the given index\n");
		} else if (c == 'S') {
			verify_strassen(&a, &cnt, &vek_n, &vek_m);
		} else if (c == 'Q') {
			deallocate(&a, &vek_n, &vek_m, cnt, &v_sum);
			return 0;
		} else if (strchr("ABEGHIJKLNQUVXYZ", c)) {
			printf("Unrecognized command\n");
		}
	}
	return 0;
}
