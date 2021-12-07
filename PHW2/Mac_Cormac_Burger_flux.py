# Author : Bo Yuan You 
# Date 2021/11/30
# CFD program homework 2 
# To solve Burger eq with Mac Cormack Scheme MUSCL

import numpy as np
import matplotlib.pyplot as plt
import time
from RK_Burger_flux import RK_Flux_Burger_eq

class Mac_Cormac_Flux_Burger(RK_Flux_Burger_eq):
    def __init__(self, CFL, M, alpha=0.75, k=0, limiter=False) -> None:
        super().__init__(CFL, M, alpha=alpha, k=k, limiter=limiter)
        self.f = np.zeros(M+2)
        self.f1 = np.zeros(M+2)
    def solver(self, t) -> None: 
        print(f'Mac Cormac solver: t = {t:.3f}')
        if self.limiter:
            for i in range(1,self.M+1):
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
                
        self.f[-1] = self.f[1]
        self.f[0] = self.f[-2] 

        for i in range(self.M):
            self.u1[i+2] = self.u0[i+2] - self.DT/self.DX*(self.f[i+2] - self.f[i+1])
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
            
        self.f1[-1] = self.f1[1]
        self.f1[0] = self.f1[-2] 

        for i in range(self.M):
            self.u2[i+2] = self.u1[i+2] - self.DT/self.DX*(self.f1[i+1] - self.f1[i])

        self.u = 0.5 * (self.u2 + self.u0).copy()

def main():
    start = time.time()
    M20 = Mac_Cormac_Flux_Burger(0.8, 40, k=1/3, limiter=False)
    M20.full_process()
    M40 = Mac_Cormac_Flux_Burger(0.8, 40, k=1/3, limiter=True)
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
    plt.axis('auto')
    

    plt.subplot(1,2,2)
    plt.title('M = 40, CFL = 0.8, with limiter')
    plt.plot(M40.Xgrid, M40.u01, 'b-o', label= 'T = 0.1')
    plt.plot(M40.Xgrid, M40.u05, 'r-o', label= 'T = 0.5')
    plt.plot(M40.Xgrid, M40.u, 'g-o', label= 'T = 1.0')
    plt.xlabel('x')
    plt.ylabel('u(x,t)')
    plt.grid(True)
    plt.legend()
    plt.axis('auto')

    plt.savefig("images/Mac_Cormac_Flux_Burger.png")  
    plt.show()

if __name__ == '__main__':
    main()
