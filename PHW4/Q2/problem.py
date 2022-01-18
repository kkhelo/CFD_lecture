# Author : Bo Yuan You 
# Date 2022/01/15
# CFD program homework 4
# 2D MUSCL Scheme 

import numpy as np
import matplotlib.pyplot as plt


class problem():
    def __init__(self, k : float, m :float):

        self.k = k
        dx = dy = 4/m
        self.dx, self.dy = dx, dy
        self.m = m

        # central grid
        self.x = np.arange(-2-1.5*dx, 2+2.5*dx, dx)
        self.y = np.arange(-2-1.5*dy, 2+2.5*dy, dy)
        self.ul = np.zeros((m, m+1))
        self.ur = np.zeros((m, m+1))
        self.vl = np.zeros((m+1, m))
        self.vr = np.zeros((m+1, m)) 
        self.F = np.zeros((m, m+1))
        self.G = np.zeros((m+1, m))

        # corner grid
        self.xh = np.arange(-2, 2+1*dx, dx)
        self.yh = np.arange(-2, 2+1*dy, dy)

        # u
        self.u = np.zeros((self.x.size, self.y.size))
        for i in range(self.xh.size):
            for j in range(self.yh.size):
                self.u[i, j] = 1 if self.xh[i]**2 + self.yh[j]**2 <= 1 else 0
        
    def call_BC(self)->None:
        self.u[0:2:1, ::] = self.u[-4:-2:1, ::]
        self.u[-2::1, ::] = self.u[2:4:1, ::]
        self.u[::, 0:2:1] = self.u[::, -4:-2:1]
        self.u[::, -2::1] = self.u[::, 2:4:1]

    def minmod(self, U1, U2, U3)->None:       
        for i, value in enumerate(U1):
            if (value > 0) and (U2[i] > 0) and (U3[i] > 0):
                U1[i] =  min(U1[i], U2[i], U3[i])
            elif (value < 0) and (U2[i] < 0) and (U3[i] < 0):
                U1[i] =  max(U1[i], U2[i], U3[i])
            else:
                U1[i] = 0.0
        
    def get_interface(self):
        k1 = 1 - self.k
        k2 = 1 + self.k
        u = self.u
        ul = 0.25 * (k1*(u[2:-2, 1:-2] - u[2:-2, 0:-3]) + k2 * (u[2:-2, 2:-1] - u[2:-2, 1:-2]))
        ur = 0.25 * (k1*(u[2:-2, 3:] - u[2:-2, 2:-1]) + k2 * (u[2:-2, 2:-1] - u[2:-2, 1:-2]))
        vr = 0.25 * (k1*(u[1:-2, 2:-2] - u[0:-3, 2:-2]) + k2*(u[2:-1, 2:-2] - u[1:-2, 2:-2]))
        vl = 0.25 * (k1*(u[3:, 2:-2] - u[2:-1, 2:-2]) + k2*(u[2:-1, 2:-2] - u[1:-2, 2:-2]))

        for i,_ in enumerate(ul):
            self.minmod(ul[i, :], u[i, 1:-2] - u[i, 0:-3], u[i, 2:-1] - u[i, 1:-2])
            self.minmod(ur[i, :], u[i, 3:] - u[i, 2:-1], u[i, 2:-1] - u[i, 1:-2])
            self.minmod(vr[:, i], u[1:-2, i] - u[0:-3, i], u[2:-1, i] - u[1:-2, i])
            self.minmod(vr[:, i], u[3:, i] - u[2:-1, i], u[2:-1, i] - u[1:-2, i]) 

        self.ul = u[2:-2, 1:-2] + ul
        self.ur = u[2:-2, 2:-1] - ur
        self.vl = u[1:-2, 2:-2] + vl
        self.vr = u[2:-1, 2:-2] - vr

    def get_flux(self):
        alpha_x = np.maximum(self.yh, 0.01)
        alpha_y = np.maximum(self.xh, 0.01)
        for i in range(self.m):
            self.F[i, :] = 0.5*(-np.multiply(self.yh, self.ul[i, :] + self.ur[i, :]) - np.multiply(alpha_x, self.ur[i, :]-self.ul[i, :]))
            self.G[:, i] = 0.5*(-np.multiply(self.xh, self.vr[:, i] + self.vl[:, i]) - np.multiply(alpha_y, self.vr[:, i] - self.vl[:, i]))

    def solver(self, T : float, dt : float):
        t = 0.00
        while t < T:
            t = round(t + dt, 8)
            print(t)
            u0 = self.u.copy()

            # first stage
            self.get_interface()
            self.get_flux()
            self.u[2:-2, 2:-2] = self.u[2:-2, 2:-2] - dt/self.dx*(self.F[:, 1:] - self.F[:, 0:-1]) - dt/self.dy*(self.G[1:, :] - self.G[0:-1, :])
            self.call_BC()

            # second stage
            self.get_interface()
            self.get_flux()
            self.u[2:-2, 2:-2] = self.u[2:-2, 2:-2] - dt/self.dx*(self.F[:, 1:] - self.F[:, 0:-1]) - dt/self.dy*(self.G[1:, :] - self.G[0:-1, :])
            self.u = self.u + u0
            self.call_BC()

        
def main():
    problem1 = problem(k=1/3, m=64)
    problem1.solver(0.2, 0.001)
    X, Y = np.meshgrid(problem1.x, problem1.y)

    print(problem1.u[10:-10, 10:-10])
    plt.figure(figsize=(12, 9))
    plt.title('phi contour')
    contour = plt.contourf(X[2:-2, 2:-2],Y[2:-2, 2:-2],problem1.u[2:-2, 2:-2])
    plt.colorbar(contour)
    plt.savefig('images/m64dt0001.png')
    plt.show()


if __name__ == '__main__':
    main()