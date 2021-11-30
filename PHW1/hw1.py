## Author : Bo Yuan You 
## Date 2021/10/01
## CFD homework 2 
## Lecture p6~p7
## Finite Difference Method on periodic B.C.

import numpy as np
import matplotlib.pyplot as plt
# Step0 : Set dimension 
m = 20
x = np.arange(-0.1, 1.1, 0.05)
u = np.zeros(m+5)
u0 = np.zeros(m+5)

# def callBC and FDM function in advance
def callBC():
    u[0] = u[m]
    u[1] = u[m+1]
    u[m+3] = u[3]
    u[m+4] = u[4]

def FDM(arr1, arr2):
    for i in range(m+1):
        arr1[i+2] = arr2[i+2] - (DT/DX*(arr2[i+2] - arr2[i+1]))

# Step1 : Grid
DX = 1/m
CFL = 0.8
DT = DX * CFL
for i in range(x.size):
    x[i] = (i-2) * DX

# Step2 : Initial data  
for i in range(m+1):
    u[i+2] = np.sin(2 * np.pi * x[i+2])

# Step3: B.C.
callBC()

# Step4 : PDE(FDM), main part of the program 
t = 0
plt.plot(x,u,'r-*', label = 'u(x,0)')
while(round(t + DT, 2) <= 2):
    t = round(t + DT, 2)    # Round off the two decimal places
    print(t)
    u0 = u.copy()
    FDM(u,u0) 
    callBC()


# Step5 : Plotting 
plt.plot(x,u,'b-o', label = 'u(x,2)')
plt.title("u(x,0) vs. u(x,2)")
plt.xlabel("x")
plt.ylabel("u(x,t)")
plt.grid(True)
plt.legend()
plt.savefig("hw1.png")
plt.show()
