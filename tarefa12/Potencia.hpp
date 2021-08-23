#include<vector>
#include <iostream>
#include<cmath>
using namespace std;

class Potencia {
  private:
    vector<double> v0;
    vector<vector<double>> A;
    vector<double> x1;
    vector<vector<double>> matrizL;
    vector<vector<double>> matrizU;
    double dominante;
    double erro;
  public:
    Potencia(vector<vector<double>> A);
    double potenciaRegular();
    double potenciaInversa();

    //Métodos Auxiliares
    void normalizar(vector<double> &vetor);
    vector<double> multiplicar(vector<vector<double>> A, vector<double> x);
    double produtoEscalar(vector<double> v1, vector<double> v2);
    vector<double> fatoracaoLU(vector<vector<double>> A, vector<double> x);
    void iniciarMatrizL();
    void iniciarMatrizU();
    vector<double> pivoteamento(vector<double> linhaPivo, vector<double> linhaOperanda, double multiplicador);
    vector<double> solucionarSistemaInferior(vector<vector<double>> matriz, vector<double> x);
    vector<double> solucionarSistemaSuperior(vector<vector<double>> matriz, vector<double> x);
};

Potencia::Potencia(vector<vector<double>> A) {
  this->v0 = {1,1,1};
  this->A = A;
  this->erro = 0.0001;
  this->dominante = 0;
  this->potenciaRegular();
}

double Potencia::potenciaRegular() {
  vector<double> vetNovo = this->v0;
  double domNovo = 0;
  do{
    this->dominante = domNovo;
    this->v0 = vetNovo;
    normalizar(vetNovo);
    this->x1 = vetNovo;
    vetNovo = multiplicar(A, vetNovo);
    domNovo = produtoEscalar(x1,vetNovo);
  }while((domNovo - this->dominante)/domNovo > this->erro);

  return domNovo;
}

double Potencia::potenciaInversa() {
  
 
 
 
  return 0;
}

void Potencia::normalizar(vector<double> &vetor) {
  double norma = 0;
  for(int i=0; i<vetor.size(); i++) {
    norma += vetor[i]*vetor[i];
  }
  norma = sqrt(norma);
  for(int i=0; i<vetor.size(); i++) {
    vetor[i] = vetor[i]*(1/norma);
  }
}

vector<double> Potencia::multiplicar(vector<vector<double>> A, vector<double> x) {
  vector<double> v;
  double temp;
  for(int i=0; i<A.size(); i++) {
    for(int j=0; j<A.size(); j++) {
      temp += A[i][j]*x[j];
    }
    v.push_back(temp);
    temp = 0;
  }
  return v;
}

double Potencia::produtoEscalar(vector<double> v1, vector<double> v2) {
  double produto = 0;
  for(int i=0; i<v1.size(); i++) {
    produto += v1[i]*v2[i];
  }
  return produto;
}



void Potencia::iniciarMatrizL() {
    this->matrizL.resize(this->A.size());
    for(int i=0; i<A.size(); i++) {
        this->matrizL[i].resize(this->A.size());
        for(int j=0; j<this->A.size(); j++) {
            if(i == j)  this->matrizL[i][j] = 1;
            else        this->matrizL[i][j] = 0;
        }
    }
}

void Potencia::iniciarMatrizU() {
  matrizU.resize(this->A.size());
  for(int i=0; i<this->A.size(); i++) {
      matrizU[i].resize(this->A.size());
      for(int j=0; j<this->A.size(); j++) {
          matrizU[i][j] = this->A[i][j];
      }
  }
}

vector<double> Potencia::pivoteamento(vector<double> linhaPivo, vector<double> linhaOperanda, double multiplicador) {
  for(int j=0;j<this->A.size(); j++) {
        linhaOperanda[j] = linhaOperanda[j] - multiplicador*linhaPivo[j];
    }
    return linhaOperanda;
}

vector<double> Potencia::solucionarSistemaInferior(vector<vector<double>> matriz, vector<double> f) {
  vector<double> x(matriz.size());
      double controlador;
      for(int i=0;i<x.size();i++) {
          controlador = 0;
          for(int j=0; j<i; j++) {
              if(matriz[i][j] != 0) {
                  controlador = controlador + matriz[i][j]*x[j];  
              }
          }
          x[i] = (f[i]-controlador)/matriz[i][i];
      }
      return x;
}

vector<double> Potencia::solucionarSistemaSuperior(vector<vector<double>> matriz, vector<double> f) {
    vector<double> x(matriz.size());
    double controlador;
    for(int i=x.size()-1; i>=0; i--) {
        controlador=0;
        for(int j=x.size()-1; j>i; j--) {
            if(matriz[i][j] != 0) {
                controlador = controlador + matriz[i][j]*x[j];  
            }
        }
        x[i] = (f[i]-controlador)/matriz[i][i];
    }
    return x;
}

vector<double> Potencia::fatoracaoLU(vector<vector<double>> A, vector<double> x) {
    int iteracao = 0;
    double pivo;
    double multiplicador;

    for(int j=0; j<A.size(); j++) {
        pivo = matrizU[j][j];
        for(int i=j+1; i<A.size(); i++) {
            multiplicador = matrizU[i][j]/pivo;
            matrizL[i][j] = multiplicador;
            matrizU[i] = this->pivoteamento(matrizU[iteracao],matrizU[i], multiplicador);
            multiplicador = matrizU[i][j]/pivo;
        }
        iteracao++;
    }
    vector<double> y(this->A.size());
    y = solucionarSistemaInferior(this->matrizL, x);
    vector<double> novo(A.size());
    novo = this->solucionarSistemaSuperior(this->matrizU, y);
    for(int i=0; i<novo.size(); i++) {
        novo[i] = abs(novo[i]);
    }
    return novo;
}