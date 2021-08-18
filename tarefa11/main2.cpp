#include<iostream>



using namespace std;

 
int main()
{
    Matrix2f m;
    m << 2,-5,
         3, 1;
   
    Vector2f v;
    v << 15,31;

    cout << (m.transpose() * m).ldlt().solve(m.transpose() * v) << endl;
    // Matrix2f y = Matrix2f::Random();
    // cout << "Here is the matrix m:" << endl << m << endl;
    // cout << "Here is the matrix y:" << endl << y << endl;
    // Matrix<float,3,2> x = m.fullPivLu().solve(y);
    // if((m*x).isApprox(y))
    // {
    // cout << "Here is a solution x to the equation mx=y:" << endl << x << endl;
    // }
    // else
    // cout << "The equation mx=y does not have any solution." << endl;

}