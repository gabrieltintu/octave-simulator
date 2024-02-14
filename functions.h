// Tintu Gabriel-Claudiu 2023 - 2024

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int mod(int n);
void allocate_rows_collumns(int **vek_n, int **vek_m, int cnt);
void allocate_list(int ****a, int n, int m, int cnt);
void allocate_rows_colls_red(int **vek_l, int **vek_c, int *nlin, int *ncol);
void print_matrices(int x, int ***a, int *vek_n, int *vek_m, int cnt);
void red_mat(int ****a, int x, int n, int m, int *vek_l, int *vek_c, int nr);
void multiply_mat(int ****a, int x, int y, int **vek_n, int **vek_m, int cnt);
void matrix_sum(int ***a, int **v_sum, int *vek_n, int *vek_m, int cnt);
void sort_matrices(int ****a, int **v_sum, int **vek_n, int **vek_m, int cnt);
void transpose_matrix(int ****a, int **vek_n, int **vek_m, int x);
int **mat_product(int **mat_aux1, int **mat_aux2, int n, int ***rez);
int **power_matrix(int ***a, int n, int p, int x, int ***rez);
void free_matrix(int ****a, int **vek_n, int **vek_m, int x, int *cnt);
void deallocate(int ****a, int **vek_n, int **vek_m, int cnt, int **v_sum);
void verify_R(int x, int p, int *vek_n, int *vek_m, int cnt, int ****a);
int **allocate_matrix(int n);
int **add(int **aux1, int **aux2, int n);
int **substract(int **aux1, int **aux2, int n);
void free_strassen(int **m1, int **m2, int **m3, int **m4, int n);
void values_1(int ***sa1, int ***sb1, int **aux1, int **aux2, int n);
void values_2(int ***sa2, int ***sb2, int **aux1, int **aux2, int n);
void values_3(int ***sa3, int ***sb3, int **aux1, int **aux2, int n);
void values_4(int ***sa4, int ***sb4, int **aux1, int **aux2, int n);
int **val_mat(int **mat, int **cm1, int **cm2, int **cm3, int **cm4, int n);
int **strassen2(int **aux1, int **aux2, int n);
void verify_strassen(int ****a, int *cnt, int **vek_n, int **vek_m);
