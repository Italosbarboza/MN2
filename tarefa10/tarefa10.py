from math import sqrt

# Size Vetor
def normOfVector(vec):
    comp = 0
    vectAux = vec
    for j in range(3):
            comp = comp + vectAux[j]*vec[j]
    comp = sqrt(comp)
    return comp


def unitVector(vec):
    sum = 0

    size = len(vec)
    
    sum = normOfVector(vec)
    vectAux = vec
    for i in range(size):
        vec[i] = vectAux[i]/sum
    return vectAux

def dotProduct(v1, v2):
    size = len(v1)
    s = 0
    for i in range(size):
        s = s + v1[i]*v2[i]
    return s


def eigenvector(matrix, vnew, vold):
    size = len(vnew)
    
    x = len(matrix)
    power = 0
    lambNew = 1
    lambOld = 0
    eps = 0.0000001
    step = 0
    
    # numberMultiplications = input('Enter the number of multiplications: ')

    # power = int(numberMultiplications)

    while(abs((lambNew - lambOld)/lambNew) > eps):
        step = step + 1
        lambOld = lambNew

        vnew = unitVector(vnew)
        vAux = vnew
        vold = vAux

        vAux2 = []
        for k in range(size):
            summ = 0
            for z in range(size):
                summ = summ + matrix[k][z]*vold[z]
            vAux2.append(summ);

        vnew = vAux2;

        lambNew = dotProduct(vold, vnew)

    return("lambda", step, "=", lambNew, "e o auto-vetor =", vold)

     
print('Para a Matriz [5, 2, 1], [2, 3, 1], [1, 1, 2]')
matrix1 = [[5, 2, 1], [2, 3, 1], [1, 1, 2]]
matrix2 = [[40, 8, 4, 2, 1], [8, 30, 12, 6, 2], [4, 12, 20, 1, 2], [2, 6, 1, 25, 4], [1, 2, 2, 5, 5]]
vnew2 = [1, 0, 0, 0, 0]
vold2 = [0, 0, 0, 0, 0]

vnew = [1, 0, 0]
vold = [0, 0, 0]
normv = 0
power = 0

a = eigenvector(matrix1, vnew, vold)
print(a)

print('Para a Matriz [40, 8, 4, 2, 1], [8, 30, 12, 6, 2], [4, 12, 20, 1, 2], [2, 6, 1, 25, 4], [1, 2, 2, 5, 5]')

b = eigenvector(matrix2, vnew2, vold2)
print(b)
