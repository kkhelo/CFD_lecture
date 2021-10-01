import numpy as np
import matplotlib.pyplot as plt


m = 20
x = np.arange(-0.1, 1.1, 0.05)
u = np.zeros(m+5)
u0 = np.zeros(m+5)

DX = 1/m
CFL = 0.8
DT = DX * CFL

def callBC():
    u[0] = u[m]
    u[1] = u[m+1]
    u[m+3] = u[3]
    u[m+4] = u[4]

def FDM(arr1, arr2):
    for i in range(m+1):
        arr1[i+2] = arr2[i+2] - (DT/DX*(arr2[i+2] - arr2[i+1]))


for i in range(x.size):
    x[i] = (i-2) * DX

for i in range(m+1):
    u[i+2] = np.sin(2 * np.pi * x[i+2])

callBC()

t = 0
while(t <= 2):
    print(t)
    u0 = u.copy()
    FDM(u,u0)
    callBC()
    if(t == 0):
        plt.plot(x,u,'r-*', label = 'u(x,0)')
    t = round(t + DT, 2)


plt.plot(x,u,'b-*', label = 'u(x,2)')
plt.title("u(x,0) compare to  u(x,2)")
plt.xlabel("x")
plt.xlabel("u(x,t)")
plt.grid(True)
plt.legend()
plt.savefig("hw1.png")
plt.show()

