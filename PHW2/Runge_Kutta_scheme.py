# Author : Bo Yuan You 
# Date 2021/11/30
# CFD program homework 2 
# First ordr upwind scheme for 1D linear convection 

import numpy as np
import matplotlib.pyplot as plt
import time
from euler_scheme import Euler_explicit_scheme

class Runge_Kutta_two_steps(Euler_explicit_scheme):
    def __init__(self, CFL, M) -> None:
        super().__init__(CFL, M)
        self.u1 = np.zeros(M+4)
        self.u2 = np.zeros(M+4)
        self.f1 = np.zeros(M+1)

    def solver(self, t) -> None:
        print(f'Using Runge-Kutta solver: t = {t:.2f}')
        for i in range(self.M+1):
            self.f[i] = 0.5 * (self.u0[i+1] + self.u0[i+2])

        for i in range(self.M):
            self.u1[i+2] = self.u0[i+2] - self.DT/self.DX*(self.f[i+1] - self.f[i])
        self.call_BC(self.u1)

        for i in range(self.M+1):
            self.f1[i] = 0.5 * (self.u1[i+1] + self.u1[i+2])
        
        for i in range(self.M):
            self.u2[i+2] = self.u1[i+2] - self.DT/self.DX*(self.f1[i+1] - self.f1[i])

        self.u = 0.5 * (self.u2 + self.u0).copy()

def main():
    start = time.time()
    M20 = Runge_Kutta_two_steps(0.8, 20)
    M20.full_process()
    M40 = Runge_Kutta_two_steps(0.8, 40)
    M40.full_process()
    end = time.time() - start 
    print('Time pass : ', end)
    plt.figure()
    plt.plot(M20.Xgrid, M20.u, 'b-o', label= 'M = 20')
    plt.plot(M40.Xgrid, M40.u, 'r-o', label= 'M = 40')
    plt.xlabel('x')  
    plt.ylabel('u(x,t)')
    plt.title('Two Step Runge Kutta Method')
    plt.grid(True)
    plt.legend()
    plt.axis([-0.1, 1.1, -1.2, 1.2])
    plt.savefig("images/Runge_Kutta_two_steps.png")  
    plt.show()

if __name__ == "__main__":
    main()