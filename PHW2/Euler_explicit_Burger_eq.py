# Author : Bo Yuan You 
# Date 2021/11/30
# CFD program homework 2 
# To solve Burger eq with Euler explict method

import numpy as np
import matplotlib.pyplot as plt
import time
from euler_scheme import Euler_explicit_scheme


class Euler_scheme_Burger_eq(Euler_explicit_scheme):
    def __init__(self, CFL, M, alpha = 0.75) -> None:
        super().__init__(CFL, M)
        self.f = np.zeros(self.M+5)
        self.alpha = alpha

    def init_data(self) -> None:
        self.u = 0.5 + 0.25*np.sin(2*np.pi*self.Xgrid)
        self.u01 = np.zeros(self.M+5)
        self.u05 = np.zeros(self.M+5)

    def solver(self, t) -> None:
        print(f'Using Euler explicit solver : t = {t:.3f}')
        self.f = (np.square(self.u0)/2).copy()
        print(self.f)
        for i in range(self.M+1):
            self.u[i+2] = self.u0[i+2] - self.DT/self.DX/2*(self.f[i+3] - self.f[i+1] - self.alpha*(self.u0[i+3] + self.u0[i+1]))
        
    def PDE_solver(self) -> None:
        t = 0
        while t < 1.0:
            t = round(t + self.DT, 3)
            self.u0 = self.u.copy()
            self.solver(t)
            self.call_BC(self.u)
            if (t == 0.1) or (t == 0.5):
                if t == 0.1:
                    self.u01 = self.u.copy()
                    print('u01 build')
                else:
                    self.u05  = self.u.copy()
                    print('u05 build')

            
def main():
    start = time.time()
    M20 = Euler_scheme_Burger_eq(0.5, 20)
    M20.full_process()
    M40 = Euler_scheme_Burger_eq(0.8, 40)
    M40.full_process()
    end = time.time() - start 
    print(M20.u, M20.u05)
    print('Time pass : ', end)

    plt.figure()
    # plt.title('Burger Equation with Euler Explicit Method')
    plt.plot(M20.Xgrid, M20.u01, 'b-o', label= 'T = 0.1')
    plt.plot(M20.Xgrid, M20.u05, 'b-o', label= 'T = 0.5')
    plt.plot(M20.Xgrid, M20.u, 'b-o', label= 'T = 1.0')
    plt.title('M = 20, CFL = 0.5')
    plt.xlabel('x')
    plt.ylabel('u(x,t)')
    plt.grid(True)
    plt.legend()
    plt.savefig("images/Euler_scheme_burger_M20.png") 
    plt.show()
    

    plt.figure()
    # plt.title('Burger Equation with Euler Explicit Method')
    plt.plot(M40.Xgrid, M40.u01, 'r-o', label= 'T = 0.1')
    plt.plot(M40.Xgrid, M40.u05, 'r-o', label= 'T = 0.5')
    plt.plot(M40.Xgrid, M40.u, 'r-o', label= 'T = 1.0')
    plt.title('M = 40, CFL = 0.8')
    plt.xlabel('x')
    plt.ylabel('u(x,t)')
    plt.grid(True)
    plt.legend()

    plt.axis([-0.1, 1.1, -1.2, 1.2])
    plt.savefig("images/Euler_scheme_burger_M40.png")  
    plt.show()

if __name__ == '__main__':
    main()