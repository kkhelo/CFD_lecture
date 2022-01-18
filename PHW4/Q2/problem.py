# Author : Bo Yuan You 
# Date 2022/01/15
# CFD program homework 4
# 2D MUSCL Scheme 

import numpy as np


class problem1():
    def __init__(self, k : float, dx : float, dy : float):

        self.k = k
        self.dx = dx
        self.dy = dy

        # central grid
        self.x = np.arange(-2-0.5*dx, 2+0.5*dx, dx)
        self.y = np.arange(-2-0.5*dy, 2+0.5*dy, dy)
        self.ul = np.zeros((self.x.size, self.y.size))
        self.ur = np.zeros((self.x.size, self.y.size))
        self.vl = np.zeros((self.x.size, self.y.size))
        self.vr = np.zeros((self.x.size, self.y.size))

        # corner grid
        self.xh = np.arange(-2-2*dx, 2+3*dx, dx)
        self.yh = np.arange(-2-2*dx, 2+3*dy, dy)

        # u
        self.u = np.zeros((self.xh.size, self.yh.size))
        for i in range(self.xh.size):
            for j in range(self.yh.size):
                self.u[i, j] = 1 if self.xh[i]**2 + self.yh[j]**2 <= 1 else 0
        
    def call_BC(self)->None:
        self.u[0:2:1, ::] = self.u[-5:-3:1, ::]
        self.u[-2::1, ::] = self.u[3:5:1, ::]
        self.u[::, 0:2:1] = self.u[::, -5:-3:1]
        self.u[::, -2::1] = self.u[::, 3:5:1]

    def minmod(self, U1, U2, U3)->None:       
        for i, value in enumerate(U1):
            if (value > 0) and (U2[i] > 0) and (U3[i] > 0):
                U1[i] =  min(U1[i], U2[i], U3[i])
            elif (value < 0) and (U2[i] < 0) and (U3[i] < 0):
                U1[i] =  max(U1[i], U2[i], U3[i])
            else:
                U1[i] = 0.0
        
    def get_interface(self, k):
        k1 = 1 - k
        k2 = 1 + k
        u = self.u
        self.ul = 0.25 * (k1*(u[1:-2:1, ::] - u[0:-3:1, ::]) + k2 * (u[2:-1:1, ::] - u[1:-2:1, ::]))
        self.ur = 0.25 * (k1*(u[3::1, ::] - u[2:-1:1, ::]) + k2 * (u[2:-1:1, ::] - u[1:-2:1, ::]))
        self.vr = 0.25 * (k1*(u[::, 1:-2:1] - u[::, 0:-3:1]) + k2*(u[::, 2:-1:1] - u[::, 1:-2:1]))
        self.vl = 0.25 * (k1*(u[::, 3::1] - u[::, 2:-1:1]) + k2*(u[::, 2:-1:1] - u[::, 1:-2:1]))



test = problem1(1/3, 1/8, 1/8)

print(test.ul.shape)
print(test.u.shape)