#include <vector>
#include <iostream>
#include<string>
#include<cmath>
using namespace std;

class Funcao {
    public:
        vector<double> coeficientes = {};
        void entrada();
        double f(double x);
        Funcao();
};

class NewtonCotes {
    private:
        Funcao funcao;
        double pontoInicial;
        double pontoFinal;
        double delta;
    public:

        NewtonCotes(Funcao funcao, double pontoInicial, double  pontoFinal);
        double integralFechadaGrau1(double inicio, double final);
        double integralFechadaGrau2(double inicio, double final);
        double integralFechadaGrau3(double inicio, double final);
        double integralFechadaGrau4(double inicio, double final);
        double integralAbertaGrau1(double inicio, double final);
        double integralAbertaGrau2(double inicio, double final);
        double integralAbertaGrau3(double inicio, double final);
        double integralAbertaGrau4(double inicio, double final);
        double particionar(double inicio, double final);
};

Funcao::Funcao() {}

void Funcao::entrada() {
    string coefString;
    string coef = " ";
    cout << "Insira os coeficientes da sua funcao polinomial: " << endl;
    getline(cin,coefString);
    for(int i=0; i < coefString.size(); i++) {
        if(coefString[i] == ' ') {
            this->coeficientes.push_back(stod(coef));
            coef = ' ';
        }
        else if(i < coefString.size()-1){
            coef += coefString[i];
        }
        else {
            coef = coef + coefString[i];
            this->coeficientes.push_back(stod(coef));
        }
    }
}

double Funcao::f(double x) {
    double resultado = 0;
    for(int i=0; i<this->coeficientes.size(); i++) {
         resultado += pow(x,coeficientes.size()-i-1) * this->coeficientes[i];
    }
    return resultado;
}

NewtonCotes::NewtonCotes(Funcao funcao, double pontoInicial, double  pontoFinal) {
    this->funcao = funcao;
    this->pontoFinal = pontoFinal;
    this->pontoInicial = pontoInicial;
    this->delta = 0.0001;
    this->particionar(this->pontoInicial,this->pontoFinal);
}

double NewtonCotes::integralFechadaGrau1(double inicio, double final) {
    double h = (this->delta)/2;
    double x1 = this->funcao.f(inicio)+this->funcao.f(final);
    return h*(x1);
}

double NewtonCotes::integralFechadaGrau2(double inicio, double final) {
    double h = (this->delta)/2.00;
    double x1,x2,x3;
    double termo;
    x1 = this->funcao.f(inicio);
    x2 = this->funcao.f(inicio+h);
    x3 = this->funcao.f(final);
    
    termo = 4.00*x2;

    return ((h/3.00)*(x1+termo+x3));
}

double NewtonCotes::integralFechadaGrau3(double inicio, double final) {
    double h = this->delta/3.00;
    double x1,x2,x3,x4;
    double segundoTermo, terceiroTermo;
    x1 = this->funcao.f(inicio);
    x2 = this->funcao.f(inicio+h);
    x3 = this->funcao.f(inicio+2.00*h);
    x4 = this->funcao.f(final);
    segundoTermo = 3.00*x2;
    terceiroTermo = 3.00*x3;
    return (3.00/8.00)*h*(x1+segundoTermo+terceiroTermo+x4);
}
double NewtonCotes::integralFechadaGrau4(double inicio, double final) {
    double h = (this->delta)/4;
    double x1,x2,x3,x4;
    double primeiroTermo, segundoTermo, terceiroTermo, quartoTermo;
    x1 = this->funcao.f(inicio+h);
    x2 = this->funcao.f(inicio+2*h);
    x3 = this->funcao.f(inicio+3*h);
    x4 = this->funcao.f(inicio+4*h);

    primeiroTermo = (106.00/45.00)*x1;
    segundoTermo = (2.00/5.00)*x2;
    terceiroTermo = (26.00/15.00)*x3;
    quartoTermo = (14.00/45.00)*x4;
    return h*(primeiroTermo - segundoTermo + terceiroTermo + quartoTermo);
}

double NewtonCotes::integralAbertaGrau1(double inicio, double final) {
    double h = this->delta/3.00;
    double x1,x2;
    x1 = this->funcao.f(inicio+h);
    x2 = this->funcao.f(inicio+2*h);
    return (this->delta/2.00)*(x1+x2);
}

double NewtonCotes::integralAbertaGrau2(double inicio, double final) {
    double h = this->delta/4.00;
    double x1, x2, x3;
    double primeiroTermo, terceiroTermo;

    x1 = this->funcao.f(inicio+h);
    x2 = this->funcao.f(inicio+(2*h));
    x3 = this->funcao.f(inicio+(3*h));
    primeiroTermo = 2*x1;
    terceiroTermo = 2*x3;
    return (h*4.00/3.00)*((primeiroTermo)-x2+(terceiroTermo));
}

double NewtonCotes::integralAbertaGrau3(double inicio, double final) {
    double h = this->delta/5.00;
    double x1,x2,x3,x4;
    double primeiroTermo, quartoTermo;
    x1 = this->funcao.f(inicio+h);
    x2 = this->funcao.f(inicio+2*h);
    x3 = this->funcao.f(inicio+3*h);
    x4 = this->funcao.f(inicio+4*h);

    primeiroTermo = 11*x1;
    quartoTermo = 11*x4;
    return (5.00/24.00)*h*(primeiroTermo+x2+x3+quartoTermo);
}

double NewtonCotes::integralAbertaGrau4(double inicio, double final) {
    double h = (this->delta)/6;
    double x1,x2,x3,x4;
    double primeiroTermo, segundoTermo, terceiroTermo, quartoTermo;
    x1 = this->funcao.f(inicio+2*h);
    x2 = this->funcao.f(inicio+3*h);
    x3 = this->funcao.f(inicio+4*h);
    x4 = this->funcao.f(inicio+5*h);

    primeiroTermo = (57.00)*x1;
    segundoTermo = (21.00)*x2;
    terceiroTermo = (9.00)*x3;
    quartoTermo = (33.00)*x4;

    return (h/10.00)*(primeiroTermo - segundoTermo - terceiroTermo + quartoTermo);
}

double NewtonCotes::particionar(double inicio, double final) {
    double grau1Aberta = 0, grau1Fechada = 0,
           grau2Aberta = 0, grau2Fechada =0,
           grau3Aberta = 0, grau3Fechada =0,
           grau4Aberta = 0, grau4Fechada =0;

    cout << "Integracao Newton-Cotes no intervalo [" 
         << inicio << "," << final <<"]" << endl;
    for(double i=inicio; i<final; i+=this->delta) {
         grau1Aberta += this->integralAbertaGrau1(i,i+this->delta);
         grau1Fechada += this->integralFechadaGrau1(i,i+this->delta);
         grau2Aberta += this->integralAbertaGrau2(i,i+this->delta);
         grau2Fechada += this->integralFechadaGrau2(i,i+this->delta);
         grau3Aberta += this->integralAbertaGrau3(i,i+this->delta);
         grau3Fechada += this->integralFechadaGrau3(i,i+this->delta);
         grau4Aberta += this->integralAbertaGrau4(i,i+this->delta);
         grau4Fechada += this->integralFechadaGrau4(i,i+this->delta);
    }
     cout << "Integral na Abordagem Aberta grau 1: "<< grau1Aberta << endl;
     cout << "Integral na Abordagem Fechada grau 1: "<< grau1Fechada << endl;
     cout << "Integral na Abordagem Aberta grau 2: "<< grau2Aberta << endl;
     cout << "Integral na Abordagem Fechada grau 2: "<< grau2Fechada << endl;
     cout << "Integral na Abordagem Aberta grau 3: "<< grau3Aberta << endl;
     cout << "Integral na Abordagem Fechada grau 3: "<< grau3Fechada << endl;
     cout << "Integral na Abordagem Aberta grau 4: "<< grau4Aberta << endl;
     cout << "Integral na Abordagem Fechada grau 4: "<< grau4Fechada << endl;
}