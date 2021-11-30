# Author : Bo Yuan You 
# Date 2021/11/30
# CFD program homework 2 
# First ordr upwind scheme for 1D linear convection 

import numpy as np
import matplotlib.pyplot as plt
import time

class Euler_explicit_scheme():
    def __init__(self, CFL, M) -> None:
        self.CFL = CFL
        self.M = M        
        self.DX = 1/M
        self.DT = round(CFL * self.DX, 3)
        self.Xgrid = np.linspace(0-1.5*self.DX, 1+0.5*self.DX, self.M+3)
        self.u = np.zeros(M+3)
        self.u0 = np.zeros(M+3)
        self.f = np.zeros(M+2)

    def init_data(self) -> None:
        for i in range(self.M+1):
            self.u[i+1] = np.sin(2*np.pi*self.Xgrid[i+1])

    def call_BC(self, arr) -> None:
        arr[0] = arr[self.M]
        arr[self.M+2] = arr[2]

    def solver(self, t) -> None:
        print(f'Using Euler explicit solver : t = {t:.2f}')
        for i in range(self.M+2):
            self.f[i] = 0.5 * (self.u[i] + self.u[i+1])
        for i in range(self.M + 1):
            self.u[i+1] = self.u0[i+1] - (self.DT/self.DX*(self.f[i+1] - self.f[i]))

    def PDE_solver(self) -> None:
        t = 0
        while t < 2.0:
            t = round(t + self.DT, 2)
            self.u0 = self.u.copy()
            self.solver(t)
            self.call_BC(self.u)

    def full_process(self) -> None:
        self.init_data()
        self.call_BC(self.u)
        self.PDE_solver()
        print('')

def main():
    start = time.time()
    M20 = Euler_explicit_scheme(0.8, 20)
    M20.full_process()
    M40 = Euler_explicit_scheme(0.8, 40)
    M40.full_process()
    end = time.time() - start 
    print('Time pass : ', end)
    plt.figure()
    plt.plot(M20.Xgrid, M20.u, 'b-o', label= 'M = 20')
    plt.plot(M40.Xgrid, M40.u, 'r-o', label= 'M = 40')
    plt.xlabel('x')
    plt.ylabel('u(x,t)')
    plt.title('Euler Explicit Scheme')
    plt.grid(True)
    plt.legend()
    # plt.axis([-0.1, 1.1, -1.2, 1.2])
    plt.savefig("images/Euler_scheme.png")  
    plt.show()



if __name__ == "__main__":
    main()
    # test = Euler_explicit_scheme(0.8, 20)
    