#include<vector>
#include<iostream>
#include<cmath>

using namespace std;

class LU {
    private:
        vector<vector<double>> matrizL;
        vector<vector<double>> matrizU;
        vector<vector<double>> A;
        vector<double> b;

    public:
        LU(vector<vector<double>> A, vector<double> b);
        vector<double> fatoracaoLU(vector<vector<double>> A, vector<double> x);
        void iniciarMatrizL();
        void iniciarMatrizU();
        vector<double> pivoteamento(vector<double> linhaPivo, vector<double> linhaOperanda, double multiplicador);
        vector<double> solucionarSistemaInferior(vector<vector<double>> matriz, vector<double> b);
        vector<double> solucionarSistemaSuperior(vector<vector<double>> matriz, vector<double> b);
};

LU::LU(vector<vector<double>> A, vector<double> b) {
    this->A = A;
    this->b = b;
    this->iniciarMatrizL();
    this->iniciarMatrizU();
}

vector<double> LU::fatoracaoLU(vector<vector<double>> A, vector<double> x) {
    LU(A,x);
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
    return novo;
}
void LU::iniciarMatrizL(){
    this->matrizL.resize(this->A.size());
    for(int i=0; i<A.size(); i++) {
        this->matrizL[i].resize(this->A.size());
        for(int j=0; j<this->A.size(); j++) {
            if(i == j)  this->matrizL[i][j] = 1;
            else        this->matrizL[i][j] = 0;
        }
    }
}
void LU::iniciarMatrizU(){
    matrizU.resize(this->A.size());
    for(int i=0; i<this->A.size(); i++) {
      matrizU[i].resize(this->A.size());
      for(int j=0; j<this->A.size(); j++) {
          matrizU[i][j] = this->A[i][j];
      }
  }
}
vector<double> LU::pivoteamento(vector<double> linhaPivo, vector<double> linhaOperanda, double multiplicador){
    for(int j=0;j<this->A.size(); j++) {
        linhaOperanda[j] = linhaOperanda[j] - multiplicador*linhaPivo[j];
    }
    return linhaOperanda;
}
vector<double> LU::solucionarSistemaInferior(vector<vector<double>> matriz, vector<double> b){
    vector<double> x(matriz.size());
    double controlador;
    for(int i=0;i<x.size();i++) {
        controlador = 0;
        for(int j=0; j<i; j++) {
            if(matriz[i][j] != 0) {
                controlador = controlador + matriz[i][j]*x[j];  
            }
        }
        x[i] = (b[i]-controlador)/matriz[i][i];
    }
    return x;
}
vector<double> LU::solucionarSistemaSuperior(vector<vector<double>> matriz, vector<double> b){
    vector<double> x(matriz.size());
    double controlador;
    for(int i=x.size()-1; i>=0; i--) {
        controlador=0;
        for(int j=x.size()-1; j>i; j--) {
            if(matriz[i][j] != 0) {
                controlador = controlador + matriz[i][j]*x[j];  
            }
        }
        x[i] = (b[i]-controlador)/matriz[i][i];
    }
    return x;
}
