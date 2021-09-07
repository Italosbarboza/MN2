import numpy as np
from functools import reduce
import math

def normalizar(vetor):
    comprimento = np.sqrt(reduce(lambda acumulador, coordenada: acumulador + coordenada**2, vetor, 0))
    return list(map(lambda coordenada : coordenada/comprimento, vetor))

def comprimento_vetor(vector):
    return np.sqrt(reduce(lambda acumulador, coordenada: acumulador + coordenada**2, vector, 0))

#Função para fazer arredondamentos nas matrizes criadas pelo NumPy
def mascara(matriz, casas_decimais):
    matriz_com_mascara = [[0] * len(matriz) for i in range (len(matriz))]
    
    for i in range(len(matriz)):
        for j in range(len(matriz)):
            matriz_com_mascara[i][j] = round(matriz[i][j], casas_decimais)
    
    return matriz_com_mascara

def matriz_de_householder_por_coluna(matriz, coluna):
    
    #Criando dois vetores nulos
    w = np.zeros(len(matriz))
    w_linha = np.zeros(len(matriz))
    
    #Copiando os elementos abaixo da diagonal da coluna i da matriz A
    for i in range (coluna, len(matriz)-1):
        w[i+1] = matriz[i+1][coluna]

    comp_w = comprimento_vetor(w)

    #Copiando comprimento de w para a posição coluna + 1 do vetor w'
    w_linha[coluna+1] = comp_w

    #Calculando e normalizando N
    N = w - w_linha
    n = normalizar(N)

    #Montando matriz de Householder
    H =  np.identity(len(matriz)) - 2*np.outer(n,n)
    
    return mascara(H,3)  

def metodo_de_householder(matriz):
    n = len(matriz)
    #Inicializando Matrizes
    H = mascara(np.identity(n, dtype=float),1)
    matriz_velha = matriz

    for i in range (0, n-2):
        #Matriz de Householder 
        H_novo = matriz_de_householder_por_coluna(matriz_velha, i)
        
        #Transformação de Similaridade
        matriz_nova = mascara(np.dot(np.dot(H_novo,matriz_velha),H_novo),2)
        
        #Salvando Matriz
        matriz_velha = matriz_nova

        #Matriz de Acumulação
        H = mascara(np.dot(H, H_novo),4)
    
    return (matriz_nova, H)

#Mascara para arredondamento
def mascara(vetor, casas_decimais):

    if len(vetor.shape) == 1:
        matriz_com_mascara = [round(item, casas_decimais) for item in vetor]

    else:        
        matriz_com_mascara = [[0] * len(vetor) for i in range (len(vetor))]
    
        for i in range(len(vetor)):
            for j in range(len(vetor)):
                matriz_com_mascara[i][j] = round(vetor[i][j], casas_decimais)
        
    return matriz_com_mascara

def matriz_jacobi_por_elemento_ij (matrix, i, j, n):
    error = 0.000001
    theta = 0.0

    matrix = np.matrix(matrix)
    J = np.identity(n)

    if abs(matrix[i,j]) <= error:
        return J

    if abs(matrix[j,j]) <= error:
        if matrix[i,j] < 0:
            theta = (math.pi)/2
        else:
            theta = - (math.pi)/2
    
    else:
        theta = math.atan( -(matrix[i,j]) / matrix[j,j])
    
    J[i,i], J[j,j], J[i,j], J[j,i] = math.cos(theta), math.cos(theta), math.sin(theta), - math.sin(theta)

    return J

def decomposicao_QR(matrix, n):
    QT = np.identity(n)
    R_velha = np.matrix(matrix)

    for j in range (n - 1):
        for i in range (j+1, n):
            J_ij = matriz_jacobi_por_elemento_ij(R_velha, i, j, n)
            R_nova = np.dot(J_ij, R_velha)
            R_velha = R_nova
            QT = np.dot(J_ij, QT)

    Q = np.transpose(QT)
    R = R_nova

    return(Q,R)

def soma_dos_quadrados_termos_abaixo_diagonal(matriz, n):
    soma = 0
    for i in range (n-1):
        for j in range(i + 1, n):
            soma += math.pow(matriz[i, j], 2)
    return soma 

def metodo_QR(matrix, n, error):
    lamb = np.zeros(n)
    vetor_acumulo_qr = []
    erro_atual = 100

    P = np.identity(n)
    matriz_velha = matrix

    while (erro_atual > error):
        (Q, R) = decomposicao_QR(matriz_velha, n)
        matriz_nova = np.dot(R,Q)
        matriz_velha = matriz_nova
        vetor_acumulo_qr.append(matriz_nova)

        P = np.dot(P, Q)
        
        erro_atual = soma_dos_quadrados_termos_abaixo_diagonal(matriz_nova, n)

    for i in range (n):
        lamb[i] = matriz_nova[i,i]
    
    return (mascara(P,6), mascara(lamb,6), vetor_acumulo_qr)

# Main

A = [[40,8,4,2,1], [8,30,12,6,2], [4,12,20,1,2],
    [2,6,1,25,4], [1,2,2,4,5]]

met_qr = metodo_QR(A, 5, 0.000001)

print('Questão 1')

print('Matriz Diagonal')
print(np.diag(met_qr[1]))


print('Matriz Acumulada P')
print(np.matrix(met_qr[0]))

'''
print('Matrizes de Cada Varredura')
for matriz in met_qr[2]:
    print(matriz, '\n')
'''

print('Pares (Autovalor, Autovetor) de A')
for i in range (5):
    autovalor = met_qr[1][i]
    autovetor = np.transpose(met_qr[0])[i]
    print('({},{})'.format(autovalor, autovetor))


print('\nQuestão 2')
(matriz_tridiagonal, matriz_householder) = metodo_de_householder(A)

#A matriz tridiagonal é passada no metodo de Jacobi
matriz_resultante = metodo_QR(matriz_tridiagonal, 5, 0.000001)[0]
array_autovalores = metodo_QR(matriz_tridiagonal, 5, 0.000001)[1]

#Para obter a matriz de autovetores correta, é preciso multiplicar a matriz de householder pela resultante do metodo de jacobi
matriz_autovetores = np.dot(matriz_householder, matriz_resultante)

print('Pares (Autovalor, Autovetor) de A do método de Householder.')
for i in range (5):
    autovalor = array_autovalores[i]
    autovetor = np.transpose(matriz_autovetores)[i]
    print('({},{})'.format(autovalor, autovetor))
