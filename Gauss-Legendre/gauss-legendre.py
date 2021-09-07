import math

# @info: com s pertencente ao conj de raizes do polinômio de Legendre de grau n.
def x(a, b, s):
    return (a + b)/2 + (b - a)/2 * s

def f(x):
    return pow(math.sin(2*x) + 4*pow(x, 2) + 3*x, 2)

def raizes_pol_legendre(n):
    if(n == 2):
        return [-0.577350, 0.577350]
    elif (n == 3):
        return [-0.774597, 0, 0.774597]
    elif (n == 4):
        return [0.861136, 0.339981, -0.339981, -0.861136]
    else:
        return None

def pesos_interpol_lagrange(n):
    if(n == 2):
        return [1, 1]
    elif (n == 3):
        return [0.555556,  0.888889, 0.555556]
    elif (n == 4):
        return [0.347854, 0.652145, 0.652145, 0.347854]
    else:
        return None

# @info: função aproximada de f utilizando o metodo de gauss-legendre
# @param:
#     a - valor inicial
#     b - valor final 
#     n - quantidade pontos de interpolação
#     f - função
#     d - quantidade de divisoes
def integral_aprox_part(n, f, d, a=0, b=1):
    s = raizes_pol_legendre(n)
    w = pesos_interpol_lagrange(n)
    delta_x = (b-a)/pow(2, d)
    
    Xi = 0
    Xf = 0
    somatorio = 0
    for j in range(int(pow(2, d))):
        Xi = a + (j * delta_x)
        Xf = Xi + delta_x

        for Sj, Wj in zip(s, w):
            somatorio += f( x(Xi, Xf, Sj) ) * Wj
    
    return (Xf-Xi)/2 * somatorio

# @info: teste de tolerância que verifica a quantidade de iterações necessárias até a função aproximada com n pontos de interpolação atingir a tolerância e
def teste_tolerancia(n, e):
    d = 1
    result_ant = integral_aprox_part(n, f, d)

    while True:
        d += 1
        result_atual = integral_aprox_part(n, f, d)
        if ( abs((result_atual - result_ant)/result_atual) < e ):
            return (d, result_atual)
        else:
            result_ant = result_atual
def main():
    e = pow(10, -6)
    
    print("### TESTE DE TOLERANCIA ###\n F: 'pow(sin(2*x) + 4*pow(x, 2) + 3*x, 2)'\n E: {}".format(e))
    for n in range(2, 5):
        d, r = teste_tolerancia(n, e)
        print("\n N: {} pontos de interpolação\n I: {} iterações\n P: {} partições\n VA~ {:.6f}".format(n, d, int(pow(2, d)), r))
    
if __name__ == '__main__':
    main()