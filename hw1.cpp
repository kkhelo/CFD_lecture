// CFD lecture program 1
#include <math.h>
#include <stdio.h>
#include <iostream>

using namespace std;

int main(void)
{

    // Initialized array size
    int m = 20;
    double x[m], u[m], u0[m];
    for (int i = -2; i <= m+2; i++)
    {
        x[i] = 0;
        u[i] = 0;
        u0[i] = 0;
    }

    // first step. grid
    double DX = 1 / double(m);
    double CFL = 0.8;
    double DT = CFL * DX;
    for (int i = -2; i <= m + 2; i++)
    {
        x[i] = i * DX;
    }

    // second step. initial data
    for (int i = 0; i <= m; i++)
    {
        u[i] = sin(2 * M_PI * x[i]);
    }

    // thrid step. B.C.
    u[-2] = u[m-2];
    u[-1] = u[m -1];
    u[m + 1] = u[1];
    u[m + 2] = u[2];

    // fourth step. PDE(FDM) Main part of the program
    for (double t = 0; t <= 2; t += DT)
    {
        // reserve known time and u will become new time
        for(int i = -2; i <= 20; i++) u0[i] = u[i];           
        // main part
        for (int i = 0; i <= m; i++)
        {
            u[i] = u0[i] - DT / DX * (u0[i] - u0[i-1]); 
        }
        // call BC again
        u[-2] = u[m-2];
        u[-1] = u[m -1];
        u[m + 1] = u[1];
        u[m + 2] = u[2];
    }

    for (int i = -2; i <= m + 2; i++)
    {
        cout << u0[i] << " ";
    }
}
