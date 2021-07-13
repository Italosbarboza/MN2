#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;

class GaussLegendre {

     public: 

          double a_, b_, tol_, N_;
          
          GaussLegendre(double a, double b, double N, double tol) {
               a_ = a;
               b_ = b;
               N_ = N;
               tol_ = tol;
               cout << "Integral de " << a << " ate " << b << "." << endl
                    << "Com " << N << " particoes e tolerancia de " << tol << "." << endl
                    << endl;
          }

          //f(x) = (sen(2x) + 4x^2+ 3x)^2
          double funcao(double x) { 
               double var = sin(2 * x) + 4 * pow(x, 2) + 3 * x;
               return pow(var, 2);
          }

          double x(double ak, double xi, double xf) {
               return (xi + xf) / 2 + ((xf - xi) / 2) * ak;
          }

          double integralTolerancia(int pontosInterpolacao, double *w, double *a) {
    
               double Iv, In = 0,
               delta, sum,
               xi, xf, tol = 1;
               int count, N = 1;

               do {
                    Iv = In;
                    In = 0;
                    count = 0;

                    delta = (b_ - a_) / N;

                    while (count < N) {
                         xi = a_ + count * delta;
                         xf = xi + delta;
                         sum = 0;

                         for (int k = 0; k < pontosInterpolacao; k++)
                              sum += (w[k] * funcao(x(a[k], xi, xf)));

                         In += ((xf - xi) / 2) * sum;
                         count++;
                    }

                    N *= 2;
                    tol = (double)fabs((In - Iv) / In);
               } while (tol > tol_);

               cout << "NUMERO DE ITERACOES: " << count << endl
               << "Integral = ";

               return In;
          }

          double integralParticao(int pontosInterpolacao, double *w, double *a) {
               double I = 0, delta = (b_ - a_) / N_,
               xi, xf, sum;

               for (int i = 0; i < N_; i++) {
                    xi = a_ + i * delta;
                    xf = xi + delta;
                    sum = 0;
                    for (int k = 0; k < pontosInterpolacao; k++)
                         sum += w[k] * funcao(x(a[k], xi, xf));
                    I += ((xf - xi) / 2) * sum;
               }

               return I;
          }

          double pontosInterpolacao2(int t) {
               double a[2], w[2];

               a[0] = -0.5773502691;
               a[1] = -a[0];

               w[0] = w[1] = 1;

               switch (t) {
                    case 1:
                         return integralParticao(2, w, a);
                    case 2:
                         return integralTolerancia(2, w, a);
               }
               return 0;
          }

          double pontosInterpolacao3(int t) {
               double a[3], w[3];

               a[0] = -0.7745966692;
               a[1] = 0;
               a[2] = -a[0];

               w[0] = w[2] = 0.5555555555;
               w[1] = 0.8888888888;

               switch (t) {
                    case 1:
                         return integralParticao(3, w, a);
                    case 2:
                         return integralTolerancia(3, w, a);
               }
               return 0;
          }

          double pontosInterpolacao4(int t) {
               double a[4], w[4];

               a[0] = 0.8611363115;
               a[1] = -a[0];
               a[2] = 0.3399810435;
               a[3] = -a[2];

               w[0] = w[1] = 0.3478548451;
               w[2] = w[3] = 0.6521451548;

               switch (t) {
                    case 1:
                         return integralParticao(4, w, a);
                    case 2:
                         return integralTolerancia(4, w, a);
               }
               return 0;
          }

};

int main() {
     cout << "Tarefa 5: Integrais Gauss-Legendre" << endl
          << "FUNCAO : f(x) = (sen(2x) + 4x^2+ 3x)^2" << endl
          << endl;

     GaussLegendre p = GaussLegendre(0, 1, 1, pow(10, -6));

     cout << p.funcao(10);

     cout << "2 PONTOS DE INTEPOLACAO: " << p.pontosInterpolacao2(1) << endl
          << "3 PONTOS DE INTEPOLACAO: " << p.pontosInterpolacao3(1) << endl
          << "4 PONTOS DE INTEPOLACAO: " << p.pontosInterpolacao4(1) << endl
     << endl;

     cout << "GAUSS LEGENDRE - TOLERANCIA: " << endl
          << endl
          << "2 PONTOS DE INTERPOLACAO" << endl
          << p.pontosInterpolacao2(2) << endl
          << endl
          << "3 PONTOS DE INTERPOLACAO" << endl
          << p.pontosInterpolacao3(2) << endl
          << endl
          << "4 PONTOS DE INTERPOLACAO" << endl
          << p.pontosInterpolacao4(2) << endl;

     return 0;
}