import numpy as np

def PVC2(n: int, fxy = 4) -> np.ndarray:
  delta_x = 1.0/n
  delta_y = 1.0/n

  mask = {
    'hborder': 1/(delta_x**2),
    'vborder': 1/(delta_y**2),
    'center': -2 * (1/(delta_x**2) + 1/(delta_y**2))
  }

  dim = (n-1)**2

  # Find indexes matrix
  index_matrix = np.empty((dim, dim), dtype=int)
  index_matrix.fill(-1)

  for i in range(dim):
    if i % (n-1) != 0: index_matrix[i][0] = i-1
    if (i+1) % (n-1) != 0: index_matrix[i][1] = i+1
    if i >= (n-1) != 0: index_matrix[i][2] = i-(n-1)
    if i < (dim-(n-1)) != 0: index_matrix[i][3] = i+(n-1)

  # Calculate coefficient matrix
  coeff_matrix = np.zeros((dim, dim), dtype=float)
  for i in range (dim):
    coeff_matrix[i][i] = mask['center']
    for j in range(4):
      if index_matrix[i][j] != -1: # get index
        k = index_matrix[i][j]
        coeff_matrix[i][k] = mask['hborder']

  # Create values vector
  values = np.empty((dim,), dtype=float)
  values.fill(fxy)
  # Calculate result
  inv_coeff_matrix = np.linalg.inv(coeff_matrix)
  result = np.dot(inv_coeff_matrix, values)
  return result

vector = PVC2(8)
tamVector = len(vector)
print("Solução Aproximada")
for i in range (tamVector):
    print('y',i+1, "=", vector[i])


