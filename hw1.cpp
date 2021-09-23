// CFD lecture program 1
#include <math.h>
#include <vector>
#include <stdio.h>
#include <iostream>

using namespace std;

int main(void)
{

    // Preallocate vector size
    int m = 20;
    vector<double> x(m + 5); // include -2, -1, 0, ..... m +1, m + 2
    vector<double> u(m + 5);
    vector<double> u0(m + 5);

    // first step. grid
    double DX = 1 / double(m);
    double CFL = 0.8;
    double DT = CFL * DX;
    for (int i = 0; i < m + 5; i++)
    {
        x[i] = (i - 2) * DX;
    }

    // second step. initial data
    for (int i = 0; i <= m; i++)
    {
        u[i + 2] = sin(2 * M_PI * x[i + 2]);
    }

    // thrid step. B.C.
    u[0] = u[m];     // u(-2) = u(m-2)
    u[1] = u[m + 1]; // u(-1) = u(m-1)
    u[m + 3] = u[3]; // u(m+1) = u(1)
    u[m + 4] = u[4]; // u(m+2) = u(2)

    // fourth step. PDE(FDM) Main part of the program
    for (double t = 0; t <= 2; t += DT)
    {
        u0 = u; // reserve known time and uwill become new time
        for (int i = 0; i <= m; i++)
        {
            u[i + 2] = u0[i + 2] - DT / DX * (u0[i + 2] - u0[i + 1]);
        }
        // call BC again
        u[0] = u[m];     // u(-2) = u(m-2)
        u[1] = u[m + 1]; // u(-1) = u(m-1)
        u[m + 3] = u[3]; // u(m+1) = u(1)
        u[m + 4] = u[4]; // u(m+2) = u(2)
    }
}
