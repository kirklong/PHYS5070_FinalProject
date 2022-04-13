#!/usr/bin/env python
from functions import *
import numpy as np

class coords:
    def __init__(self,x,y,z):
        self.x = x
        self.y = y
        self.z = z
    
    def __repr__(self):
        return "x = {0:.2f}, y = {1:.2f}, z = {2:.2f}".format(self.x,self.y,self.z)

class velocities:
    def __init__(self,vx,vy,vz):
        self.x = vx
        self.y = vy
        self.z = vz

    def __repr__(self):
        return "vx = {0:.2f}, vy = {1:.2f}, vz = {2:.2f}".format(self.x,self.y,self.z)


class body:
    def __init__(self,coords,velocities,mass):
        self.pos = coords
        self.v = velocities
        self.m = mass

    def __repr__(self):
        posStr = self.pos.__repr__()
        vStr = self.v.__repr__()
        return "{0}\n{1}\nm = {2:.2f}".format(posStr,vStr,self.m)

class system:
    def __init__(self,bodies):
        self.nBodies = len(bodies)
        self.bodies = bodies

    def update(self,Δt,G=6.67408313131313e-11):
        x,y,z,vx,vy,vz = nBodyStep(self.bodies,self.nBodies,Δt,G) #from functions.py
        for n in range(self.nBodies):
            self.bodies[n] = body([x[n],y[n],z[n]],[vx[n],vy[n],vz[n]],self.bodies[n].m)






