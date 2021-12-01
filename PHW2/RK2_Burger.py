# Author : Bo Yuan You 
# Date 2021/11/30
# CFD program homework 2 
# To solve Burger eq with Runge Kutta scheme

import numpy as np
import matplotlib.pyplot as plt
import time
from Euler_explicit_Burger_eq import Euler_scheme_Burger_eq

class RK_Burger(Euler_scheme_Burger_eq):
    def __init__(self, CFL, M, alpha = 0.75) -> None:
        super().__init__(CFL, M)
        self.alpha = alpha
        self.u1 = np.zeros(M+5)
        self.u2 = np.zeros(M+5)
        self.f1 = np.zeros(M+4)

    def solver(self, t) -> None:
        print(f'Using Runge-Kutta solver: t = {t:.2f}')
        for i in range(self.M+4):
            self.f[i] = 0.5*(self.u[i]**2/2 + self.u[i+1]**2/2 - self.alpha*(self.u[i+1] - self.u[i]))

        for i in range(self.M + 1):
            self.u1[i+2] = self.u0[i+2] - self.DT/self.DX*(self.f[i+2] - self.f[i+1])
        self.call_BC(self.u1)

        for i in range(self.M+4):
            self.f1[i] = 0.5*(self.u1[i]**2/2 + self.u1[i+1]**2/2 - self.alpha*(self.u1[i+1] - self.u1[i]))
        
        for i in range(self.M + 1):
            self.u2[i+2] = self.u1[i+2] - self.DT/self.DX*(self.f1[i+2] - self.f1[i+1])

        self.u = 0.5 * (self.u2 + self.u).copy()

def main():
    start = time.time()
    M20 = RK_Burger(0.5, 20)
    M20.full_process()
    M40 = RK_Burger(0.8, 40)
    M40.full_process()
    end = time.time() - start 
    print('Time pass : ', end)

    plt.figure(figsize=(12,5))
    plt.title('Burger Equation with Runge Kutta Method')

    plt.subplot(1,2,1)
    plt.title('M = 20, CFL = 0.5')
    plt.plot(M20.Xgrid, M20.u01, 'b-o', label= 'T = 0.1')
    plt.plot(M20.Xgrid, M20.u05, 'r-o', label= 'T = 0.5')
    plt.plot(M20.Xgrid, M20.u, 'g-o', label= 'T = 1.0')
    plt.xlabel('x')
    plt.ylabel('u(x,t)')
    plt.grid(True)
    plt.legend()
    plt.axis([0, 1, 0, 1])
    

    plt.subplot(1,2,2)
    plt.title('M = 40, CFL = 0.8')
    plt.plot(M40.Xgrid, M40.u01, 'b-o', label= 'T = 0.1')
    plt.plot(M40.Xgrid, M40.u05, 'r-o', label= 'T = 0.5')
    plt.plot(M40.Xgrid, M40.u, 'g-o', label= 'T = 1.0')
    plt.xlabel('x')
    plt.ylabel('u(x,t)')
    plt.grid(True)
    plt.legend()
    plt.axis([0, 1, 0, 1])

    plt.savefig("images/Runge_Kutta_burger.png")  
    plt.show()

if __name__ == '__main__':
    main()