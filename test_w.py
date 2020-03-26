# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 22:13:37 2020

@author: Juliano
"""
import numpy as np

class DiffusingField:#i think inherit Solver?
    def __init__(self, 
                 D = None, gamma = None,
                 xdim = None, ydim = None, time_iters = None,
                 dx = None, dy = None, dt = None,
                 BC_type = None, bc_x = None, bc_y = None):
        
        self.current_step = 0
        
        self._init_grid(dt,dx,dy,time_iters,xdim,ydim)
        
        self._init_diff_params(D,gamma,BC_type,bc_x,bc_y)
                        
        self.concentrations = self._init_concentration_field(
                                len(self.xdim),len(self.ydim),self.itrs_t)
        
    def _init_diff_params(self,D,gamma,BC,bc_x,bc_y):
        self.D = D if D is not None else .1*self.dx*self.dx/self.dt# (D*dt/dx**2) << 1
        self.g = gamma if gamma is not None else 0
        self.gamma = self.g
        
        if BC is None:
            self.BC = 'floating'
        else:
            self.BC = BC
        
        
        
        # self._set_select_BC(BC,bc_x,bc_y)
        
    # def _set_select_BC(self,BC,bc_x,bc_y):
    def _init_grid(self,dt,dx,dy,time_iters,xdim,ydim):
        dt = dt if dt is not None else .00005
        dx = dx if dx is not None else .01
        dy = dy if dy is not None else .01
        
        self.dx = dx
        self.dy = dy
        self.dt = dt
        
        if xdim is not None:
            self.xdim = np.arange(0,xdim,dx)
        else:
            self.xdim = np.arange(0,.1,dx)
        if ydim is not None:
            self.ydim = np.arange(0,ydim,dx)
        else:
            self.ydim = np.array([0])
        
        if time_iters is not None:
            self.time = np.arange(0,time_iters,dt)
            self.itrs_t = time_iters
        else:
            self.time = np.arange(0,10,dt)
            self.itrs_t = 11
        
        return
    
    def _init_concentration_field(self,x,y,t):
        if y != 1:
            field = np.zeros(shape = (x,y,t))
            
            return field
        else:
            field = np.zeros(shape = (x,t))
            self.step = self._step_1d
            return field
        return field
    
    # def _calc_next_step_1d(self,next_, curr, plus, minus, z, alpha, beta):
    #     next_ = curr + alpha * (plus + minus) - beta * z
    #     return next_
    def _step_1d(self):
        
        # dt = dt if dt is not None else .00005
        # dx = dx if dx is not None else .01
        # D =.49*self.dx*self.dx/self.dt
        curr_step = self.current_step
        D = self.D
        g = self.g
        
        dt = self.dt
        dx = self.dx
        
        c =np.zeros(self.concentrations.shape)
        c[:,:] = self.concentrations[:,:]
        
        for ix in range(1,len(c)-1):
            c[ix,curr_step+1] = (c[ix,curr_step] + (D*dt/dx**2) *(
            c[ix+1,curr_step] + c[ix-1,curr_step] - 2*c[ix,curr_step])
             - dt*g*c[ix,curr_step] )
        
        c[0,curr_step+1] = (c[0,curr_step] + (D*dt/dx**2) *(
            c[1,curr_step] - c[0,curr_step])
             - dt*g*c[0,curr_step] )
        c[-1,curr_step+1] = (c[-1,curr_step] + (D*dt/dx**2) *(
            c[-2,curr_step] - c[-1,curr_step])
             - dt*g*c[-1,curr_step] )
        
        
        
        self.concentrations[:,:] = c[:]
        self.current_step +=1

    def evolve(self):
        
        for itr in range(self.itrs_t-1):
            self.step()



a = DiffusingField()
