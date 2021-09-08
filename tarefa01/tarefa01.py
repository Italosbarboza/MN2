from math import sqrt, pow, e

def functionCalculate(x):
    return sqrt(pow(e, 3*x) + (4*pow(x,2)))

def secondCentralDerivativeError4(delta, x):
    return (1/pow(delta,2))*((-1/12)*functionCalculate(x+(2*delta)) + ((4/3)*functionCalculate(x + delta)) + ((4/3) * functionCalculate(x - delta)) - ((1/12) * functionCalculate(x - 2*delta)) - ((5/2) * functionCalculate(x)))

print(secondCentralDerivativeError4(0.3125, 2))
