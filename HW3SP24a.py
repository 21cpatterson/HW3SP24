import math
import copy
import random

def is_symmetric_positive_definite(matrix):
    """
    Check if the matrix is symmetric and positive definite.

    Parameters:
    - matrix: 2D list, the matrix to be checked.

    Returns:
    - bool: True if the matrix is symmetric and positive definite, False otherwise.
    """
    n = len(matrix)

    # Check if the matrix is square
    if len(matrix) != len(matrix[0]):
        return False

    # Check if the matrix is symmetric
    for i in range(n):
        for j in range(i + 1, n):
            if matrix[i][j] != matrix[j][i]:
                return False

    # Check if the matrix is positive definite
    try:
        cholesky_decomposition(matrix)
        return True
    except ValueError:
        return False

def cholesky_decomposition(matrix):
    """
    Perform Cholesky decomposition on a symmetric positive definite matrix.

    Parameters:
    - matrix: 2D list, the matrix to be decomposed.

    Returns:
    - L: 2D list, lower triangular matrix such that matrix = L * L^T.
    """
    n = len(matrix)
    L = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1):
            sum_val = sum(L[i][k] * L[j][k] for k in range(j))
            if i == j:
                L[i][j] = math.sqrt(matrix[i][i] - sum_val)
            else:
                L[i][j] = (1.0 / L[j][j]) * (matrix[i][j] - sum_val)

    return L

def forward_substitution(L, b):
    """
    Perform forward substitution to solve the system Ly = b.

    Parameters:
    - L: 2D list, lower triangular matrix.
    - b: list, the right-hand side vector.

    Returns:
    - y: list, solution vector for Ly = b.
    """
    n = len(L)
    y = [0.0] * n

    for i in range(n):
        y[i] = b[i]
        for j in range(i):
            y[i] -= L[i][j] * y[j]
        y[i] /= L[i][i]

    return y

def backward_substitution(L_transpose, y):
    """
    Perform backward substitution to solve the system L^T x = y.

    Parameters:
    - L_transpose: 2D list, transpose of the lower triangular matrix L.
    - y: list, the right-hand side vector.

    Returns:
    - x: list, solution vector for L^T x = y.
    """
    n = len(L_transpose)
    x = [0.0] * n

    for i in range(n - 1, -1, -1):
        x[i] = y[i]
        for j in range(i + 1, n):
            x[i] -= L_transpose[j][i] * x[j]
        x[i] /= L_transpose[i][i]

    return x

def doolittle_lu_decomposition(matrix):
    """
    Perform Doolittle LU decomposition on a matrix.

    Parameters:
    - matrix: 2D list, the matrix to be decomposed.

    Returns:
    - P: 2D list, permutation matrix.
    - L: 2D list, lower triangular matrix.
    - U: 2D list, upper triangular matrix.
    """
    n = len(matrix)
    P = [[float(i == j) for j in range(n)] for i in range(n)]
    L = [[0.0] * n for _ in range(n)]
    U = copy.deepcopy(matrix)

    for k in range(n - 1):
        pivot_row = max(range(k, n), key=lambda i: abs(U[i][k]))
        U[k], U[pivot_row] = U[pivot_row], U[k]
        P[k], P[pivot_row] = P[pivot_row], P[k]

        for i in range(k + 1, n):
            factor = U[i][k] / U[k][k]
            L[i][k] = factor
            for j in range(k, n):
                U[i][j] -= factor * U[k][j]

    return P, L, U

def lu_solve(P, L, U, b):
    """
    Solve the system of linear equations Ax = b using LU decomposition.

    Parameters:
    - P: 2D list, permutation matrix.
    - L: 2D list, lower triangular matrix.
    - U: 2D list, upper triangular matrix.
    - b: list, the right-hand side vector.

    Returns:
    - x: list, solution vector for Ax = b.
    """
    n = len(b)
    Pb = [0.0] * n

    for i in range(n):
        for j in range(n):
            Pb[i] += P[i][j] * b[j]

    # Solve Ly = Pb using forward substitution
    y = forward_substitution(L, Pb)

    # Solve Ux = y using backward substitution
    x = backward_substitution(U, y)

    return x

def print_solution(problem_num, solution, method):
    """
    Print the solution vector for a given problem.

    Parameters:
    - problem_num: int, the problem number.
    - solution: list, the solution vector.
    - method: str, the numerical method used.
    """
    print(f"\nProblem {problem_num} solution using {method} method:")
    print(solution)

# Problem 1
A1 = [[1, -1, 3, 2],
      [-1, 5, -5, -2],
      [3, -5, 19, 3],
      [2, -2, 3, 21]]

b1 = [15, -35, 94, 21]

if is_symmetric_positive_definite(A1):
    L1 = cholesky_decomposition(A1)
    y1 = forward_substitution(L1, b1)
    solution1 = backward_substitution(L1, y1)
    method1 = "Cholesky"
else:
    P1, L1, U1 = doolittle_lu_decomposition(A1)
    solution1 = lu_solve(P1, L1, U1, b1)
    method1 = "Doolittle"

print_solution(1, solution1, method1)

# Problem 2
A2 = [[4, 2, 4, 0],
      [2, 2, 3, 2],
      [4, 3, 6, 3],
      [0, 2, 3, 9]]

b2 = [20, 36, 60, 122]

if is_symmetric_positive_definite(A2):
    L2 = cholesky_decomposition(A2)
    y2 = forward_substitution(L2, b2)
    solution2 = backward_substitution(L2, y2)
    method2 = "Cholesky"
else:
    P2, L2, U2 = doolittle_lu_decomposition(A2)
    solution2 = lu_solve(P2, L2, U2, b2)
    method2 = "Doolittle"

print_solution(2, solution2, method2)
