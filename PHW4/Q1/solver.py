# Author : Bo Yuan You 
# Date 2022/01/15
# CFD program homework 4
# sover : RK MUSCL, Mac-CorMack, Jameson

import numpy as np

class RK():
    def __init__(self, U, CFL)->None:
        self.U = U
        self.CFL = CFL
        self.U_temp = np.copy(U)    
        self.eslon = 0.000001 

    def first_stage(self, F)->None:
        self.U_temp[::, 2:-3:1] = self.U[::, 2:-3:1] - self.CFL * (F[::, 1:-1:1] - F[::, 0:-2:1])
        self.U_temp[::, 0] = self.U_temp[::, 2]
        self.U_temp[::, 1] = self.U_temp[::, 2]
        self.U_temp[::, -1] = self.U_temp[::, -3]
        self.U_temp[::, -2] = self.U_temp[::, -3]
        for i in [0, 2]:
            for j,_ in enumerate(self.U_temp[i]):
                self.U_temp[i, j] = max(self.U_temp[i, j], self.eslon)
        

    def second_stage(self, F)->None:
        self.U_temp[::, 2:-3:1] = self.U_temp[::, 2:-3:1] - self.CFL * (F[::, 1:-1:1] - F[::, 0:-2:1])
        self.U_temp[::, 0] = self.U_temp[::, 2]
        self.U_temp[::, 1] = self.U_temp[::, 2]
        self.U_temp[::, -1] = self.U_temp[::, -3]
        self.U_temp[::, -2] = self.U_temp[::, -3]
        for i in [0, 2]:
            for j,_ in enumerate(self.U_temp[i]):
                self.U_temp[i, j] = max(self.U_temp[i, j], self.eslon)
        self.U = 0.5*(self.U + self.U_temp)

class MC(RK):
    def __init__(self, U, CFL) -> None:
        super().__init__(U, CFL)

    def first_stage(self, F) -> None:
        self.U_temp[::, 2:-3:1] = self.U[::, 2:-3:1] - self.CFL * (F[::, 3:-2:1] - F[::, 2:-3:1])
        self.U_temp[::, 0] = self.U_temp[::, 2]
        self.U_temp[::, 1] = self.U_temp[::, 2]
        self.U_temp[::, -1] = self.U_temp[::, -4]
        self.U_temp[::, -2] = self.U_temp[::, -4]
        self.U_temp[::, -3] = self.U_temp[::, -4]
        
        for i in [0, 2]:
            for j,_ in enumerate(self.U_temp[i]):
                self.U_temp[i, j] = max(self.U_temp[i, j], self.eslon)

    def second_stage(self, F) -> None:
        self.U_temp[::, 2:-3:1] = self.U_temp[::, 2:-3:1] - self.CFL * (F[::, 2:-3:1] - F[::, 1:-4:1])
        self.U_temp[::, 0] = self.U_temp[::, 2]
        self.U_temp[::, 1] = self.U_temp[::, 2]
        self.U_temp[::, -1] = self.U_temp[::, -4]
        self.U_temp[::, -2] = self.U_temp[::, -4]
        self.U_temp[::, -3] = self.U_temp[::, -4]
        self.U = 0.5*(self.U + self.U_temp)