//#include "Potencia.hpp"
#include "Potencia.hpp"


int main() {

  VectorXd v0(3); 
  v0 << 1,1,1;
  double erro = 0.00001;
  double M11 = 2.2;
  double M12 = 2.2;
  double M13 = 2.2;
  double M21 = 1;
  double M22 = 1;
  double M23 = 1;
  double M31 = 1;
  double M32 = 1;
  double M33 = 1;

  MatrixXd matriz0(3,3);
  matriz0 << 2,1,0,
             2,5,3,
             0,1,6;

  MatrixXd matriz1(3,3);
  matriz1 << 5,2,1,
             2,3,1,
             1,1,2;
                              
  MatrixXd matriz2(3,3);
  matriz2 << -14,1,-2,
              1, -1, 1,
            -2, 1,-11;

  MatrixXd matriz3(5,5);
  matriz3 << 40,8,4,2,1,
             8,30,12,6,2,
             4,12,20,1,2,
             2,6,1,25,4,
             1,2,2,4,5;

  Potencia *p0 = new Potencia(matriz0, v0, erro);
  Potencia *p1 = new Potencia(matriz1, v0, erro);
  Potencia *p2 = new Potencia(matriz2, v0, erro);
  Potencia *p3 = new Potencia(matriz3, v0, erro);

  cout << "autovalor: " << p1->potenciaInversa() << endl;


  // cout << "------- Para matriz 1 -------- " << endl;
  // cout << " O autovalor mais proximo de " << M11 << "é " << p1->potenciaComDeslocamento(matriz1,v0, erro, M11) << endl;
  // cout << " O autovetor correspondente é [ ";
  // p1->printarAutovetor();
  // cout << "]" << endl << endl << endl;

  // cout << " O autovalor mais proximo de " << M12 << "é " << p1->potenciaComDeslocamento(matriz1,v0, erro, M12) << endl;
  // cout << " O autovetor correspondente é [ ";
  // p1->printarAutovetor();
  // cout << "]" << endl << endl << endl;

  // cout << " O autovalor mais proximo de " << M13 << "é " << p1->potenciaComDeslocamento(matriz1,v0, erro, M13) << endl;
  // cout << " O autovetor correspondente é [ ";
  // p1->printarAutovetor();
  // cout << "]" << endl << endl << endl;

  // cout << "------- Para matriz 2 -------- " << endl;
  // cout << " O autovalor mais proximo de " << M21 << "é " << p2->potenciaComDeslocamento(matriz1,v0, erro, M21) << endl;
  // cout << " O autovetor correspondente é [ ";
  // p2->printarAutovetor();
  // cout << "]" << endl << endl << endl;

  // cout << " O autovalor mais proximo de " << M22 << "é " << p2->potenciaComDeslocamento(matriz1,v0, erro, M22) << endl;
  // cout << " O autovetor correspondente é [ ";
  // p2->printarAutovetor();
  // cout << "]" << endl << endl << endl;

  // cout << " O autovalor mais proximo de " << M23 << "é " << p2->potenciaComDeslocamento(matriz1,v0, erro, M23) << endl;
  // cout << " O autovetor correspondente é [ ";
  // p2->printarAutovetor();
  // cout << "]" << endl << endl << endl;

  // cout << "------- Para matriz 3 -------- " << endl;
  // cout << " O autovalor mais proximo de " << M31 << "é " << p3->potenciaComDeslocamento(matriz1,v0, erro, M31) << endl;
  // cout << " O autovetor correspondente é [ ";
  // p3->printarAutovetor();
  // cout << "]" << endl << endl << endl;

  // cout << " O autovalor mais proximo de " << M32 << "é " << p3->potenciaComDeslocamento(matriz1,v0, erro, M32) << endl;
  // cout << " O autovetor correspondente é [ ";
  // p3->printarAutovetor();
  // cout << "]" << endl << endl << endl;

  // cout << " O autovalor mais proximo de " << M33 << "é " << p3->potenciaComDeslocamento(matriz1,v0, erro, M33) << endl;
  // cout << " O autovetor correspondente é [ ";
  // p3->printarAutovetor();
  // cout << "]" << endl << endl << endl;

  return 0;
}