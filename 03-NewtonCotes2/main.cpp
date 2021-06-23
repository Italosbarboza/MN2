#include"newtonCotes.hpp"
int main() {
    double pontoInicial, pontoFinal;
    Funcao *funcao = new Funcao();
    funcao->entrada();
    cout << "Insira o ponto inicial do intervalo de integração: ";
    cin >> pontoInicial;
    cout << endl << "Insira o ponto final do intervalo de integração: ";
    cin >> pontoFinal;
    cout << endl;
    NewtonCotes *nc = new NewtonCotes(*funcao,pontoInicial,pontoFinal);

    return 0;
}