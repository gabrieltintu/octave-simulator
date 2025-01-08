### © 2022 Gabriel-Claudiu TINTU

# Octave Simulator

The Octave Simulator project is a C program that emulates the functionalities of Octave, a powerful programming language primarily used for numerical computations. This simulator provides a command-line interface allowing users to perform various matrix operations with ease.

## Features

- **Load Matrix (L):**
  - Allocate memory for a new matrix and load it into memory. Dynamically allocate vectors for the number of rows and columns.

- **Display Matrix Dimensions (D):**
  - Display the number of rows and columns for a specified matrix or show an error message if the matrix does not exist.

- **Print Matrix (P):**
  - Display the content of a specified matrix or show an error message if the matrix does not exist.

- **Resize Matrix (C):**
  - Resize a matrix, save it into an auxiliary matrix, free space from the matrix list, and update with the resized matrix.

- **Matrix Multiplication (M):**
  - Multiply matrices, updating the matrix list. Check for matrix existence and multiplication feasibility.

- **Order Matrices (O):**
  - Calculate the sum of each matrix, and sort matrices, row, and column vectors in ascending order based on the sum.

- **Transpose Matrix (T):**
  - Check if a matrix exists, transpose it, free the original matrix, and update with the transposed one.

- **Matrix Power (R):**
  - Check if the matrix exists, the power is positive, and the matrix is square. If conditions are met, calculate the matrix raised to the given power using a recursive function.

- **Free Matrix (F):**
  - Free the memory of a matrix with the given index and shift all matrices after it to the left by one position.

- **Strassen Matrix Multiplication (S):**
  - Check if matrices exist and can be multiplied. If conditions are met, use the Strassen algorithm to recursively calculate the multiplication, updating the matrix list.

- **Quit (Q):**
  - Deallocate all remaining allocated resources, including the matrix list, vectors for rows, columns, and sums.
