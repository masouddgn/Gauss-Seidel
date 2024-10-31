This project introduces a custom class, Data.h, which streamlines matrix-like operations on vectors.
Within this class, matrix elements are mapped to vector indices using the formula (x * Nx + y), where (x, y) represents the matrix coordinates, and Nx and Ny define the matrix dimensions.
To facilitate calculations on the grid, two interval values, hx and hy, are computed based on boundary divisions across Nx and Ny grid points.
In the primary class, we instantiate vectors u, r, and b to model the linear system represented by the equation "A.U = B".
The solution approach leverages the Gauss-Seidel iterative method, optimized with OpenMP parallelization using the Red-Black technique to enhance computational efficiency.
