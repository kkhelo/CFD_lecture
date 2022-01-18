# Author : Bo Yuan You 
# Date 2022/01/15
# CFD program homework 4
# 1D Shock problem 

import numpy as np

class problem_1D_shock():
    def __init__(self, Vl : tuple, Vr : tuple, k : float)->None:

        self.k = k

        V = np.array([np.full(52, Vl[0], dtype=float), 
                      np.full(52, Vl[1], dtype=float), 
                      np.full(52, Vl[2], dtype=float)])
        V = np.append(V, np.array([
                           np.full(53, Vr[0], dtype=float), 
                           np.full(53, Vr[1], dtype=float), 
                           np.full(53, Vr[2], dtype=float)]), axis=1)
                           
        U = np.array(V)
        U[0,::] = V[0,::]
        U[1,::] = np.multiply(V[0,::], V[1,::])
        U[2,::] = V[2,::]/(1.4-1) + 0.5 * np.divide(U[1,::]**2, U[0, ::])
        self.U = U
        self.eslon = 0.000001
        
        
    def get_U_interface(self, U :np.array)->np.array:   # get U_i+1 for MUSCL
        
        k1 = 1 - self.k
        k2 = 1 + self.k 
        Ul = np.zeros((U.shape[0], U.shape[1]-3)) 
        Ur = np.zeros((U.shape[0], U.shape[1]-3))
        SOl = np.zeros((U.shape[0], U.shape[1]-3))
        SOr = np.zeros((U.shape[0], U.shape[1]-3))


        SOl = k2 * U[::, 2:-1:1] + (k1 - k2) * U[::, 1:-2:1] - k1 * U[::, 0:-3:1]
        SOr = k1 * U[::, 3::1] + (k2 - k1) * U[::, 2:-1:1] - k2 * U[::, 1:-2:1]
        SOl *= 0.25
        SOr *= 0.25

        for i,_ in enumerate(SOl):
            self.minmod(SOl[i, ::], U[i, 2:-1:1] - U[i, 1:-2:1], U[i, 1:-2:1] - U[i, 0:-3:1])

        for i,_ in enumerate(SOr):
            self.minmod(SOr[i, ::], U[i, 2:-1:1] - U[i, 1:-2:1], U[i, 3::1] - U[i, 2:-1:1])
    
        Ul = U[::, 1:-2:1] + SOl
        Ur = U[::, 2:-1:1] - SOr
        
        for i in [0, 2]:
            for j,_ in enumerate(Ul[i]):
                Ul[i,j] = max(Ul[i,j], self.eslon)
                Ur[i,j] = max(Ur[i,j], self.eslon)

        return Ul, Ur


    def compute_flux(self, Ul : np.array, Ur : np.array)->np.array:     #get F_i+1 for MUSCL
        
        F = np.zeros(Ul.shape)
        alpha = np.zeros(Ul.shape[1])
        alphal = np.zeros(alpha.shape)
        alphar = np.zeros(alpha.shape)
        Pl = np.zeros(Ul.shape[1])
        Pr = np.zeros(Ul.shape[1])
        ul = np.zeros(Ul.shape[1])
        ur = np.zeros(Ul.shape[1])
        Hl = np.zeros(Ul.shape[1])
        Hr = np.zeros(Ul.shape[1])

        Pl = (1.4 - 1) * (Ul[2, ::] - 0.5 * (np.divide(Ul[1, ::]**2, Ul[0, ::])))
        Pr = (1.4 - 1) * (Ur[2, ::] - 0.5 * (np.divide(Ur[1, ::]**2, Ur[0, ::])))

        for i,_ in enumerate(Pl):
            Pl[i] = max(Pl[i], self.eslon)
            Pr[i] = max(Pr[i], self.eslon)

        ul = np.divide(Ul[1, ::], Ul[0, ::])
        ur = np.divide(Ur[1, ::], Ur[0, ::])
        Hl = np.divide(Ul[2, ::] + Pl, Ul[0, ::])
        Hr = np.divide(Ur[2, ::] + Pr, Ur[0, ::])
        alphal = np.sqrt(1.4*np.divide(Pl, Ul[0, ::])) + np.abs(ul)
        alphar = np.sqrt(1.4*np.divide(Pr, Ur[0, ::])) + np.abs(ur)
        alpha = np.maximum(alphal, alphar)
        

        F[0, ::] = 0.5 * (Ul[1, ::] + Ur[1, ::] - np.multiply(alpha, Ur[0, ::] - Ul[0, ::]))                     
        F[1, ::] = 0.5 * (np.multiply(ul, Ul[1, ::]) + Pl + np.multiply(ur, Ur[1, ::]) + Pr
                                - np.multiply(alpha, Ur[1, ::] - Ul[1, ::]))
        F[2, ::] = 0.5 * (np.multiply(Ul[1, ::], Hl) + np.multiply(Ur[1, ::], Hr) - np.multiply(alpha, Ur[2, ::] - Ul[2, ::]))

        return F

        
    def minmod(self, U1, U2, U3):       
        for i, value in enumerate(U1):
            if (value > 0) and (U2[i] > 0) and (U3[i] > 0):
                U1[i] =  min(U1[i], U2[i], U3[i])
            elif (value < 0) and (U2[i] < 0) and (U3[i] < 0):
                U1[i] =  max(U1[i], U2[i], U3[i])
            else:
                U1[i] = 0.0

    
    def MC_get_F(self, U : np.array)->np.array:
        F = np.zeros(U.shape, dtype=float)
        F[0, ::] = U[1, ::]
        F[1, ::] = np.divide(F[0,::]**2, U[0, ::])
        P = (1.4 - 1) * (U[2, ::] - 0.5 * F[1, ::])
        F[1, ::] = F[1, ::] + P
        H = np.divide(U[2, ::] + P, U[0, ::])
        F[2, ::] = np.multiply(F[0, ::], H)

        return F
                 