from math import sqrt

# Comprimento do Vetor
def normOfVector(vec):
    comp = 0
    vectAux = vec
    for j in range(3):
            comp = comp + vectAux[j]*vec[j]
    comp = sqrt(comp)
    return comp


def unitVector(vec):
    sum = 0
    
    sum = normOfVector(vec)
    vectAux = vec
    for i in range(3):
        vec[i] = vectAux[i]/sum
    return vectAux

def dotProduct(v1, v2):
    return(v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])

def matrixXvector(matrix, vector):
    s = 0
    for z in range(3):
        s = s + matrix[z]*vector[z]
    return s


def eigenvector(matrix, vnew, vold):

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
        #print(vnew)
        vAux = vnew
        vold = vAux


        for k in range(3):
            vnew[k] = matrixXvector(matrix[k], vold)

        print(vnew)

        #lambNew = dotProduct(vold, vnew)

        #print("lambda", step, "=", lambNew)

        #for j in range(3):
         #   vnew[j] = vnew[j]/vnew1
        #print(vnew)

matrix = [[5, 2, 1], [2, 3, 1], [1, 1, 2]]
vnew = [1, 0, 0]
vold = [0, 0, 0]
normv = 0
power = 0

eigenvector(matrix, vnew, vold)

