# Author : Bo Yuan You 
# Date 2022/01/19
# CFD program homework 5
# 2D cavity flow    
#   1) Mac CorMack scheme


from operator import mod
from tkinter.messagebox import NO
import numpy as np
import matplotlib.pyplot as plt


class MC():         # Mac CorMack
    def __init__(self, beta, dx, dy, dt, shape, nu : float = 0) -> None:
        self.beta = 1
        self.dx = dx
        self.dy = dy
        self.dt = dt
        self.nu = nu
        self.E1 = np.zeros(shape)
        self.E2 = np.zeros(shape)
        self.E1 = np.zeros(shape)
        self.F1 = np.zeros(shape)
        self.F2 = np.zeros(shape)
        self.F3 = np.zeros(shape)

    def first_stage(self, u, v, P)->None:

        E1 = u**2 + P
        E2 = np.multiply(u, v)
        E3 = self.beta * u
        F1 = np.multiply(u, v)
        F2 = v**2 + P
        F3 = self.beta * v

        a = (u[1:-1, 2:] - 2*u[1:-1, 1:-1] + u[1:-1, 0:-2]) / (self.dx**2)
        b = (u[0:-2, 1:-1] - 2*u[1:-1, 1:-1] + u[2:, 1:-1]) / (self.dy**2)
        c = (v[1:-1, 2:] - 2*v[1:-1, 1:-1] + v[1:-1, 0:-2]) / (self.dx**2)
        d = (v[0:-2, 1:-1] - 2*v[1:-1, 1:-1] + v[2:, 1:-1]) / (self.dy**2)
        u[1:-1, 1:-1] = u[1:-1, 1:-1] - self.dt/self.dx*(E1[1:-1, 2:] - E1[1:-1, 1:-1]) - self.dt/self.dy*(F1[0:-2, 1:-1] - F1[1:-1, 1:-1]) + self.dt*self.nu*(a + b)
        v[1:-1, 1:-1] = v[1:-1, 1:-1] - self.dt/self.dx*(E2[1:-1, 2:] - E2[1:-1, 1:-1]) - self.dt/self.dy*(F2[0:-2, 1:-1] - F2[1:-1, 1:-1]) + self.dt*self.nu*(c + d)
        P[1:-1, 1:-1] = P[1:-1, 1:-1] - self.dt/self.dx*(E3[1:-1, 2:] - E3[1:-1, 1:-1]) - self.dt/self.dy*(F3[0:-2, 1:-1] - F3[1:-1, 1:-1])


    def second_stage(self, u, v, P)->None:

        E1 = u**2 + P
        E2 = np.multiply(u, v)
        E3 = self.beta * u
        F1 = np.multiply(u, v)
        F2 = v**2 + P
        F3 = self.beta * v

        a = (u[1:-1, 2:] - 2*u[1:-1, 1:-1] + u[1:-1, 0:-2]) / (self.dx**2)
        b = (u[0:-2, 1:-1] - 2*u[1:-1, 1:-1] + u[2:, 1:-1]) / (self.dy**2)
        c = (v[1:-1, 2:] - 2*v[1:-1, 1:-1] + v[1:-1, 0:-2]) / (self.dx**2)
        d = (v[0:-2, 1:-1] - 2*v[1:-1, 1:-1] + v[2:, 1:-1]) / (self.dy**2)
        u[1:-1, 1:-1] = u[1:-1, 1:-1] - self.dt/self.dx*(E1[1:-1, 1:-1] - E1[1:-1, 0:-2]) - self.dt/self.dy*(F1[1:-1, 1:-1] - F1[2:, 1:-1]) + self.dt*self.nu*(a + b)
        v[1:-1, 1:-1] = v[1:-1, 1:-1] - self.dt/self.dx*(E2[1:-1, 1:-1] - E2[1:-1, 0:-2]) - self.dt/self.dy*(F2[1:-1, 1:-1] - F2[2:, 1:-1]) + self.dt*self.nu*(c + d)
        P[1:-1, 1:-1] = P[1:-1, 1:-1] - self.dt/self.dx*(E3[1:-1, 1:-1] - E3[1:-1, 0:-2]) - self.dt/self.dy*(F3[1:-1, 1:-1] - F3[2:, 1:-1])


class cavity():
    def __init__(self, m : float, Re : float, beta : float):
        self.nu = 1/Re
        self.m = m
        self.beta = beta
        dx = dy = 1/m
        x = np.arange(0 - 1*dx, 1 + 2*dx, dx)
        y = np.arange(1 + 1*dy, 0 - 2*dy, -dy)
        self.X, self.Y = np.meshgrid(x, y)
        self.dx = dx
        self.dy = dy

        # Initial condition
        self.u = np.zeros(self.X.shape)
        self.v = np.zeros(self.X.shape)
        self.P = np.ones(self.X.shape)


    def call_BC(self):

        # left
        self.u[:, 0] = -self.u[:, 2]
        self.v[:, 0] = -self.v[:, 2]
        self.P[:, 0] = self.P[:, 2]
        self.P[:, 1] = self.P[:, 2]
        self.u[:, 1] = self.v[:,1] = 0
        
        # right 
        self.u[:, -1] = -self.u[:, -3]
        self.v[:, -1] = -self.v[:, -3]
        self.P[:, -1] = self.P[:, -3]
        self.P[:, -2] = self.P[:, -3]
        self.u[:, -2] = self.v[:, -2] = 0

        #bottom 
        self.u[-1, :] = -self.u[-3, :]
        self.v[-1, :] = -self.v[-3, :]
        self.P[-1, :] = self.P[-3, :]
        self.P[-2,:] = self.P[-3, :]
        self.u[-2, :] = self.v[-2, :] = 0

        #top
        self.u[0:2, :] = self.P[0:2, :] = 1 
        
    def two_stage_solver(self, solver : MC)->None:
        
        u0 = self.u.copy()
        v0 = self.v.copy()
        P0 = self.P.copy()

        solver.first_stage(self.u, self.v, self.P)
        self.call_BC()
        solver.second_stage(self.u, self.v, self.P)
        self.call_BC()
        self.u = 0.5 * (self.u + u0)
        self.v = 0.5 * (self.v + v0)
        self.P = 0.5 * (self.P + P0)
        self.call_BC()



def main():
    model = cavity(20, 100, 10)
    model.call_BC()

    t = 0.0
    T = 100
    dt = 0.001
    solver = MC(model.beta, model.dx, model.dy, dt, model.u.shape, model.nu)

    while(t < T):
        t = round(t + dt, 5)
        print(t)

        model.two_stage_solver(solver)
    
    print(model.P)
    
    np.save('P.npy', model.P)
    np.save('u.npy', model.u)
    np.save('v.npy', model.v)

   
    

if __name__ == '__main__':
    main()