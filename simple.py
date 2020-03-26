# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 11:11:35 2020

@author: Juliano
"""
#C(i,j+1) = C(i,j) + (dt/dx^2)*(C(i+1,j) - 2*C(i,j) + C(i-1,j));
import numpy as np

def step(c, D, g, dx, dt):
    nc = np.zeros(c.shape)
    for idx in range(1,len(c)-1):
        nc[idx] = c[idx] + (D*dt/dx**2) * ( c[idx+1] + c[idx-1] - 2*c[idx] ) - dt*g*c[idx]
    
    
    nc[0] = c[0] + (D*dt/dx**2) * ( c[1] - c[0] )- dt*g*c[idx]
    nc[-1] = c[-1] + (D*dt/dx**2) * (  c[-2] - c[-1] )- dt*g*c[idx]
    return nc
def evolve(c,D,g,dx,dt):
    
    for idt in range(c.shape[1]-1):
        c[:,idt+1] = step(c[:,idt],D,g,dx,dt)
        if np.any( c < 0 ):
            
            raise Exception("negative number!")
            break
            return
    

xmax = 1
dx = 1/100
x = np.arange(0,xmax,dx)
iter_t =100
dt = 0.00005
# t = np.linspace()
v = 1/2 * dx**2/dt
D = 1. # D < 1/2 dx**2/dt

if D > v:
    raise Exception("invalid D, D must <= {}".format( v))


r = dt
rD = D*dt/dx**2 # < 1/2



gamma = 0#.00000000000000000000000000000001


c = np.zeros(shape = (len(x),iter_t))
c[5,0] = 100

evolve(c,D,gamma,dx,dt)

for col in c.T:
    print(col.sum())
    


