# Author : Bo Yuan You 
# Date 2021/11/30
# CFD program homework 2 
# To solve Burger eq with two step RK method (flux)


import numpy as np
import matplotlib.pyplot as plt
import time
from RK2_Burger import RK_Burger

class RK_Flux_Burger_eq(RK_Burger):
    def __init__(self, CFL, M, alpha=0.75, k=0, limiter=False) -> None:
        super().__init__(CFL, M, alpha=alpha)
        self.kl = 1-k
        self.kr = 1+k
        self.limiter=limiter

    def minmoid(self,a,b,c)->float:
        if a>0 and b>0 and c>0:
            return min(a,b,c)
        elif a<0 and b<0 and c<0:
            return max(a,b,c)
        else:
            return 0.0

    def solver(self, t) -> None:
        print(f'Using Runge-Kutta solver: t = {t:.2f}')

        if self.limiter:
            for i in range(self.M+1):
                du_left = self.u0[i+1] - self.u0[i]
                du_mid = self.u0[i+2] - self.u0[i+1]
                du_right = self.u0[i+3] - self.u0[i+2]
                so_left = 0.25*(self.kl*(self.u0[i+1] - self.u0[i]) + self.kr*(self.u0[i+2] - self.u0[i+1]))
                so_right = 0.25*(self.kl*(self.u0[i+2] - self.u0[i+1]) + self.kr*(self.u0[i+3] - self.u0[i+2]))
                ul = self.u0[i+1] + self.minmoid(du_left, du_mid, so_left)
                ur = self.u0[i+2] + self.minmoid(du_mid, du_right, so_right)
                self.f[i] = 0.5*(ul**2/2 + ur**2/2 - self.alpha*(ur - ul))

        else:
            for i in range(self.M+1):
                ul = self.u0[i+1] + 0.25*(self.kl*(self.u0[i+1] - self.u0[i]) + self.kr*(self.u0[i+2] - self.u0[i+1]))
                ur = self.u0[i+2] + 0.25*(self.kl*(self.u0[i+2] - self.u0[i+1]) + self.kr*(self.u0[i+3] - self.u0[i+2]))
                self.f[i] = 0.5*(ul**2/2 + ur**2/2 - self.alpha*(ur - ul))
                

        for i in range(self.M):
            self.u1[i+2] = self.u0[i+2] - self.DT/self.DX*(self.f[i+1] - self.f[i])
        self.call_BC(self.u1)

        if self.limiter:
            for i in range(self.M+1):
                du_left = self.u1[i+1] - self.u1[i]
                du_mid = self.u1[i+2] - self.u1[i+1]
                du_right = self.u1[i+3] - self.u1[i+2]
                so_left = 0.25*(self.kl*(self.u1[i+1] - self.u1[i]) + self.kr*(self.u1[i+2] - self.u1[i+1]))
                so_right = 0.25*(self.kl*(self.u1[i+2] - self.u1[i+1]) + self.kr*(self.u1[i+3] - self.u1[i+2]))
                ul = self.u1[i+1] + self.minmoid(du_left, du_mid, so_left)
                ur = self.u1[i+2] + self.minmoid(du_mid, du_right, so_right)
                self.f1[i] = 0.5*(ul**2/2 + ur**2/2 - self.alpha*(ur - ul))
        else:
            for i in range(self.M+1):
                ul = self.u1[i+1] + 0.25*(self.kl*(self.u1[i+1] - self.u1[i]) + self.kr*(self.u1[i+2] - self.u1[i+1]))
                ur = self.u1[i+2] + 0.25*(self.kl*(self.u1[i+2] - self.u1[i+1]) + self.kr*(self.u1[i+3] - self.u1[i+2]))
                self.f1[i] = 0.5*(ul**2/2 + ur**2/2 - self.alpha*(ur - ul))
            
        for i in range(self.M):
            self.u2[i+2] = self.u1[i+2] - self.DT/self.DX*(self.f1[i+1] - self.f1[i])

        self.u = 0.5 * (self.u2 + self.u0).copy()

    def PDE_solver(self) -> None:
        self.u01 = np.zeros(self.M+4)
        self.u05 = np.zeros(self.M+4)
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
    M20 = RK_Flux_Burger_eq(0.8, 40, k=1/3, limiter=False)
    M20.full_process()
    M40 = RK_Flux_Burger_eq(0.8, 40, k=1/3, limiter=True)
    M40.full_process()
    end = time.time() - start 
    print('Time pass : ', end)

    plt.figure(figsize=(12,5))
    plt.title('Burger Equation with Runge Kutta Method')

    plt.subplot(1,2,1)
    plt.title('M = 40, CFL = 0.8')
    plt.plot(M20.Xgrid, M20.u01, 'b-o', label= 'T = 0.1')
    plt.plot(M20.Xgrid, M20.u05, 'r-o', label= 'T = 0.5')
    plt.plot(M20.Xgrid, M20.u, 'g-o', label= 'T = 1.0')
    plt.xlabel('x')
    plt.ylabel('u(x,t)')
    plt.grid(True)
    plt.legend()
    plt.axis([0, 1, 0.2, 0.8])
    

    plt.subplot(1,2,2)
    plt.title('M = 40, CFL = 0.8, with limiter')
    plt.plot(M40.Xgrid, M40.u01, 'b-o', label= 'T = 0.1')
    plt.plot(M40.Xgrid, M40.u05, 'r-o', label= 'T = 0.5')
    plt.plot(M40.Xgrid, M40.u, 'g-o', label= 'T = 1.0')
    plt.xlabel('x')
    plt.ylabel('u(x,t)')
    plt.grid(True)
    plt.legend()
    plt.axis([0, 1, 0.2, 0.8])

    plt.savefig("images/Runge_Kutta_burger_flux.png")  
    plt.show()

if __name__ == '__main__':
    main()

    