#include "LU.hpp"
#include "eigen-3.3.7/Eigen/Dense"

using namespace std;
using namespace Eigen;

class Potencia {
  private:
    VectorXd v0;
    MatrixXd matriz;
    double dominante;
    VectorXd autovetor;
    double erro;
  public:
    Potencia(MatrixXd matriz, VectorXd v0, double e);
    double potenciaRegular();
    double potenciaRegular(MatrixXd matriz, VectorXd v0, double e);
    double potenciaInversa();
    double potenciaInversa(MatrixXd matriz, VectorXd v0, double e);
    double potenciaComDeslocamento(MatrixXd &matriz, VectorXd v0, double erro, double M);
    VectorXd getAutovetor();
    //Métodos Auxiliares.
    void printarVetor(VectorXd vetor);
    void normalizar(VectorXd &vetor);
    VectorXd multiplicar(MatrixXd A, VectorXd x);
    double produtoEscalar(VectorXd v1, VectorXd v2);
    void calcularTransposta(VectorXd &matriz);
    double getMaior(VectorXd x);
};



Potencia::Potencia(MatrixXd matriz, VectorXd v0, double erro) {
  this->v0 = v0;
  this->matriz = matriz;
  this->erro = erro;
  this->dominante = 0;
}

double Potencia::potenciaRegular() {
  VectorXd vetorVelho = this->v0;
  VectorXd vetorNovo = this->v0;
  double domVelho = 0;
  double domNovo;
  do {
    domVelho = domNovo;
    vetorVelho = vetorNovo;
    vetorNovo = this->multiplicar(this->matriz,vetorVelho);
    domNovo = this->produtoEscalar(multiplicar(this->matriz,vetorNovo),vetorNovo)/
              this->produtoEscalar(vetorNovo,vetorNovo);

  }while(abs((domNovo - domVelho)/domNovo) > this->erro);
  double temp = vetorNovo.size()-1;
  for(int i=0; i<=temp; i++) {
    vetorNovo[i] /= vetorNovo[temp];
  }
  
  this->autovetor = vetorNovo;
  return domNovo;
}

double Potencia::potenciaRegular(MatrixXd matriz, VectorXd v0, double e) {
  Potencia(matriz,v0,e);
  this->potenciaRegular();
}

double Potencia::potenciaInversa() {
  double alpha;
  double valorVelho, valorNovo=0;
  VectorXd z(v0.size());
  VectorXd y = this->v0;
  VectorXd aux(y.size());

  do {
    valorVelho = valorNovo;
    cout << "y: ";
    this->printarVetor(y);
    z = (this->matriz.transpose() * this->matriz).ldlt().solve(this->matriz.transpose() * y);
    cout << "z: ";
    this->printarVetor(z);
    alpha = this->getMaior(z);
    for(int i=0; i<z.size(); i++) {
      aux[i] = z[i]/y[i];
    }
    cout << "aux ";
    this->printarVetor(aux);
    valorNovo = this->getMaior(aux);
    cout << valorNovo << endl;
    for(int i=0; i<y.size(); i++) {
      y[i] = (1.00/alpha)*z[i];
    }
    cout << "--------------------" << endl;

  }while(abs((valorNovo-valorVelho)/valorNovo)>this->erro);

  double temp = y.size()-1;
  for(int i=0; i<=temp; i++) {
    y[i] /= y[temp];
  }

  this->dominante = 1.00/valorNovo;
  this->autovetor = y;
  cout << this->autovetor << endl;
  return this->dominante;
}

double Potencia::potenciaInversa(MatrixXd matriz, VectorXd v0, double e) {
  Potencia(matriz, v0, e);
  return this->potenciaInversa();
}

double Potencia::potenciaComDeslocamento(MatrixXd &matriz, VectorXd v0, double erro, double M) {
  MatrixXd A(matriz.size(), matriz.size());
  for(int i=0; i<A.size(); i++) {
    for(int j=0; j<A.size(); j++){
      if(i==j) {
        A(i,j) = matriz(i,j) - M;
      }
      else {
        A(i,j) = matriz(i,j);
      }
    }
  }
  this->dominante = this->potenciaInversa(A, v0, erro) + M;
  return this->dominante;
}

    
void Potencia::printarVetor(VectorXd vetor) {
  for(int i=0; i<vetor.size(); i++) {
    cout << vetor[i] << "  ";
  }
  cout << endl;
}

void Potencia::normalizar(VectorXd &vetor) {
  double norma = 0;
  for(int i=0; i<vetor.size(); i++) {
    norma += vetor[i]*vetor[i];
  }
  norma = sqrt(norma);
  for(int i=0; i<vetor.size(); i++) {
    vetor[i] = vetor[i]*(1/norma);
  }
}

VectorXd Potencia::multiplicar(MatrixXd A, VectorXd x) {
  VectorXd v(A.size(),x.size());
  double temp = 0;
  for(int i=0; i<A.size(); i++) {
    for(int j=0; j<A.size(); j++) {
      temp += A(i,j)*x[j];
    }
    v(i) = temp;
    temp = 0;
  }
  return v;
}

double Potencia::produtoEscalar(VectorXd v1, VectorXd v2) {
  double produto = 0;
  for(int i=0; i<v1.size(); i++) {
    produto += v1[i]*v2[i];
  }
  return produto;
}

VectorXd Potencia::getAutovetor() {
  return this->autovetor;
}

void Potencia::calcularTransposta(VectorXd &matriz) {
  MatrixXd temp = matriz;
  for(int i=0; i<matriz.size(); i++) {
    for(int j=0; j<matriz.size(); j++) {
      matriz(i,j) = temp(j,i);
    }
  }
}

double Potencia::getMaior(VectorXd x) {
  double temp=x[0];
  for(int i=0; i<x.size(); i++) {
    if(x[i]>temp)
      temp = x[i];
  }
  return temp;
}