import math
import numpy as np

#default matrix
A = [[2, 1, -1, 0],[1, 3, 0, 1], [-1, 0, 2, 1], [0, 1, 1, 4]]

#default vector
b = [1, -3, -2, -5]

#identity matrix
E = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]

#chekcing whether matrix is symmetrical
def check_def_mtrx(A):
  for i in range(4):
    for j in range(4):
      if (A[i][j] != A[j][i]):
        return False
  return True

def det_check(A):
  n = len(A)
  for i in range(n):
    if i == 1:
      if(A[0][0] < 0):
        return False
    elif i == 2:
      det = A[0][0] * A[1][1] - A[0][1] * A[1][0]
      if(det < 0):
        return False
    else:
      det = 0
      for c in range(i):
        det += ((-1) ** c) * A[0][c] * det_check(minor(A, 0, c))
      if(det < 0):
        return False   
  return True
    
def minor(matrix, row, col):
    return [row[:col] + row[col+1:] for row in (matrix[:row] + matrix[row+1:])]

#filling S and D matrix according to the formulas
def initialise_D_n_S_mtrx(S, D, A):
  #formula's base calculations same for D and S
  for i in range(4):
    for j in range(4):
      sum = 0
      for p in range(0, i):
        sum += abs(S[p][i] ** 2) * D[p][p]
      value = A[i][i] - sum

      if(i == j):
        #filling S matrix
        S[i][i] = math.sqrt(abs(value))
        #filling D matrix
        D[i][i] = np.sign(value)

      #finding S_ij elements
      elif(i < j):
        sum = 0
        for p in range(i):
          sum += S[p][i] * D[p][p] * S[p][j]
        value = (A[i][j] - sum) / ( D[i][i]*S[i][i])
        S[i][j] = value

#method call
def square_sol():
  if(check_def_mtrx(A)):
    S = np.zeros((4, 4), dtype=float)
    D = np.zeros((4, 4), dtype=float)
    initialise_D_n_S_mtrx(S, D, A)
    S_t = transposed_matrix(S)
    R = multiply_matrix(S_t, D)
    y = solve_eq_mtrx(R, b)
    x = solve_eq_mtrx(S, y)
    det = sqr_r_deteminant(D, S)
    print_results(S, D, S_t, R, y, x, det)
    find_inverse_mtrx(R, S_t)
  else:
    print("matrix is not symertrical")
    return

def Zeidel_method_sol():
  if(det_check(A)):
    print()

def print_results(S, D, S_t, R, y, x, det):
  print_matrix(S, "S Matrix")
  print_matrix(D, "D Matrix")
  print_matrix(S_t, "S^T Matrix")
  print_matrix(R, "S^T * D Matrix")
  print_vector(y, "y vector")
  print_vector(x, "x vector")
  print(f"\ndet: {det:7.2f}\n")

def find_inverse_mtrx(R, S_t):
  Y = solve_inverse_mtrx(E, R)
  X = solve_inverse_mtrx(Y, S_t)
  A_inv = transposed_matrix(X)
  print_matrix(A_inv, "A^(-1) Matrix")

def solve_inverse_mtrx(E, R):
  # L = R & U = S
  y = np.zeros(4)
  Y = np.zeros((4, 4), dtype=float)
  for j in range(4):
    y = solve_eq_mtrx(R, E[j])
    for i in range(4):
      Y[i][j] = y[i]
  return Y    

def solve_eq_mtrx(R, b):
    n = len(b)
    y = np.zeros(n)
    for i in range(n):
        sum = 0
        for j in range(i):
            sum += R[i][j] * y[j]
        y[i] = (b[i] - sum) / R[i][i]
    return y

def transposed_matrix(S):
  for i in range(4):
    for j in range(4):
        if i > j or i < i:
          temp = S[i][j]
          S[i][j] = S[j][i]
          S[j][i] = temp
  return S

def multiply_matrix(S, D):
  R = np.zeros((4, 4), dtype=float)
  for k in range(4):
    for i in range(4):
      for j in range(4):
        R[k][i] += S[k][j] * D[j][i]
  return R

def sqr_r_deteminant(D, S):
  D_multiply, S_multiply = 1, 1
  for i in range(4):
    D_multiply = D_multiply * D[i][i]
    S_multiply = S_multiply * (S[i][i] ** 2)
  determinant = D_multiply * S_multiply
  return determinant

def print_matrix(matrix, label="Matrix"):
    print(f"{label}:")
    for row in matrix:
        formatted_row = " ".join(f"{val:7.2f}" for val in row)
        print(formatted_row)
    print()

def print_vector(y, label="Vector y"):
    print(f"{label}:")
    formatted_vector = " ".join(f"{val:7.2f}" for val in y)  # Format each element in the vector
    print(formatted_vector)

if __name__ == "__main__":
  square_sol()