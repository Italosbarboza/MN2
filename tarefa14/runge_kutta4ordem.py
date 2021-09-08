import math
'''
Resumão:
    0) Fase de Inicialização
        Obter o estado S_0 pelo método de passo simples de Runge-Kutta de 2a ordem

    1) Fase de predição 
        Estimar o estado S_(i+1), i=[0, ... , n]
        Valores para Caixa Preta: S_(i-1), Si, DELTA_T, F( Si, ti ) 
            S_(i+1) = Si + DELTA_T/2  * (-1 * F(i-1) + 3 * Fi ), com
                F(i-1) = F( S_(i-1), t(i-1) )
                Fi     = F( Si    , ti )
    2) Fase de Correção (pode ser iterariva)
            S_(i+1) = Si + DELTA_T/2 * ( Fi + * F(i+1) )
            
            Para testar a convergência até a tolerância e
                S_(i+1)^(k-1) = S'(i+1)
                do 
                    ~ k++
                    S_(i+1)^k = Si + DELTA_T/2 * ( F(Si, ti) + F( S_(i+1)^(k-1), t(i+1) ) )
                    E <- ( S_(i+1)^k - S_(i+1)^(k-1) ) / S_(i+1)^k
                while (E < e) 
'''
g = 10

# F( S(t), t )
def F(k, m, v_t):
    return ( -g - (k/m)*v_t , v_t)

def range_kutta_4ordem(t_0, v_0, y_0, DELTA_T, k, m):
    i = 1
    
    t_max = t_0
    y_max = y_0

    t = t_0
    v = v_0 # S_0
    y = y_0

    # ESTIMATIVA DOS ESTADOS
    F_1 = F(k, m, v)
    # S_1 = S_0 + DELTA_T/2 * F_1
    S_1 = v + DELTA_T/2 * F_1[0]

    F_2 = F(k, m, S_1)
    # S_2 = S_0 + DELTA_T/2 * F_2
    S_2 = v + DELTA_T/2 * F_2[0]

    F_3 = F(k, m, S_2)
    # S_3 = S_0 + DELTA_T * F_3
    S_3 = v + DELTA_T * F_3[0]

    F_4 = F(k, m, S_3)
    S_4 = v + (DELTA_T/24) * (-F_1[0]*9 + F_2[0]*37 - F_3[0]*59 + F_4[0]*55)

    F_5 = F(k, m, S_4)

    # ATUALIZAÇÃO MELHORADA DO ESTADO v = Si, com i a iteração atual
    v += DELTA_T/6 * (F_2[0] + 2*F_3[0] + 2*F_3[0]+ F_4[0])
    y += DELTA_T/6 * (F_2[1] + 2*F_3[1] + 2*F_3[1]+ F_4[1])
    t += DELTA_T

    if (y > y_max):
        y_max = y
        t_max = t

    while (y > 0):
        i += 1
        
        F_1 = F_2
        F_2 = F_3
        F_3 = F_4
        F_4 = F_5  

        S_4 = v + (DELTA_T/24) * (-F_1[0]*9 + F_2[0]*37 - F_3[0]*59 + F_4[0]*55)
        F_5 = F(k, m, S_4)

        # ATUALIZAÇÃO MELHORADA DO ESTADO v = Si, com i a iteração atual
        v += DELTA_T/6 * (F_2[0] + 2*F_3[0] + 2*F_3[0]+ F_4[0])
        y += DELTA_T/6 * (F_2[1] + 2*F_3[1] + 2*F_3[1]+ F_4[1])
        t += DELTA_T

        if (y > y_max):
            y_max = y
            t_max = t

    return (y_max, t_max, t, v, i)
    

def main():
    # 1Q) Valores para PVI-2
    t_0 = 0 # s
    v_0 = 5 # m/s
    y_0 = 200 # m
    k = 0.25 # kg/s
    m = 2 # kg

    # Valores de DT da Tab. 1
    delta_ti = [0.1, 0.01, 0.001, 0.0001]
    print('-> Método de Runge-Kutta')
    for dt in delta_ti:
        res = range_kutta_4ordem(t_0, v_0, y_0, dt, k, m)
        print('\n Delta: {} seg - {} iterações\n a) altura max {}m\n b) tempo em a) {}s\n c) tempo no momento da queda ao mar {}s\n d) velocidade no momento do impacto {}m/s'.format( dt, res[4], \
                res[0], \
                res[1], \
                res[2], \
                res[3]
        ))

if __name__ == '__main__':
    main()