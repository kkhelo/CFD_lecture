# Author : Bo Yuan You 
# Date 2022/01/15
# CFD program homework 4
# 1D Shock main program

from solver import RK, MC
from problem import problem_1D_shock
import numpy as np


def main(CFL:float, T:float, dx:float)->None:
    # problem = problem_1D_shock(Vl=(1, 0, 1), Vr=(0.125, 0, 0.1), k=1/3)
    problem = problem_1D_shock(Vl=(0.44, 0.698, 3.528), Vr=(0.5, 0.0, 0.571), k=1/3)

    # initialization. compute U, F at time = 0
    U = problem.U
    
    # init solver
    # solver = RK(U, CFL=CFL)
    solver = MC(U, CFL=CFL)
    
    # start calculation 
    t = 0.0

    while t <= T:
        # print('*'*75)
        print(t)      
        t = round(t + CFL*dx, 8)

        # Ul, Ur = problem.get_U_interface(solver.U)
        # F = problem.compute_flux(Ul, Ur)
        F = problem.MC_get_F(U)

        solver.first_stage(F)

        # Ul, Ur = problem.get_U_interface(solver.U_temp)
        # F = problem.compute_flux(Ul, Ur)
        F = problem.MC_get_F(U)

        solver.second_stage(F)
        
        
    print(solver.U)

    # plotting

    import matplotlib.pyplot as plt
    plt.figure(figsize=(14,4))

    U = solver.U[::, 2:-3:1]
    x_grid = np.arange(-4.9, 5.1, 0.1)
    density = U[0, ::]
    velocity = np.divide(U[1, ::], U[0, ::])
    pressure = (1.4 - 1) * (U[2, ::] - 0.5 * (np.divide(U[1, ::]**2, U[0, ::])))

    plt.subplot(1,3,1)
    plt.plot(x_grid, density)
    plt.title('density')
    plt.subplot(1,3,2)
    plt.plot(x_grid, velocity)
    plt.title('velocity')
    plt.subplot(1,3,3)
    plt.plot(x_grid, pressure)
    plt.title('pressure')

    plt.tight_layout()
    plt.show()
    # plt.savefig('2a_Mac_CFL_0.2.png')


if __name__ == '__main__':
    main(CFL=0.5, T=2.0, dx=0.1)