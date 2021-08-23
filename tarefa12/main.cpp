#include "Potencia.hpp"

int main() {
  

  vector<vector<double>> matriz1{{5,2,1},
                                 {2,3,1},
                                 {1,1,2}};

  vector<vector<double>> matriz2 {{40,8,4,2,1},{8,30,12,6,2},{4,12,20,1,2},{2,6,1,25,4},{1,2,2,4,5}};

  Potencia *p1 = new Potencia(matriz1);
  Potencia *p2 = new Potencia(matriz2);

  cout << "O autovalor dominante para a matriz1 é " << p1->potenciaRegular() << endl;
  cout << "O autovalor dominante para a matriz2 é " << p2->potenciaRegular() << endl;
  return 0;
}