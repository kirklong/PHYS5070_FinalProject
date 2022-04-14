#!/usr/bin/env python
from functions import *
import numpy as np

class coords:
    """
    position records in cartesian space, goes like coords.x, coords.y, coords.z
    """
    def __init__(self,x,y,z):
        if isinstance(x,float) and isinstance(y,float) and isinstance(z,float):
            self.x = x
            self.y = y
            self.z = z
        else:
            print(type(x),type(y),type(z))
            raise Exception("x,y,z need to be floats")
    
    def __repr__(self):
        return "x = {0:.2f}, y = {1:.2f}, z = {2:.2f}".format(self.x,self.y,self.z)

class velocities:
    """
    velocity records in cartesian spaces, goes like velocities.x, velocities.y, velocities.z
    """
    def __init__(self,vx,vy,vz):
        if isinstance(vx,float) and isinstance(vy,float) and isinstance(vz,float):
            self.x = vx
            self.y = vy
            self.z = vz
        else:
            print(type(vx),type(vy),type(vz))
            raise Exception("vx,vy,vz need to be floats")

    def __repr__(self):
        return "vx = {0:.2f}, vy = {1:.2f}, vz = {2:.2f}".format(self.x,self.y,self.z)


class body:
    """
    class containing coords, velocities, and mass of each body (i.e. Sun, planet, etc.)
    """
    def __init__(self,p,v,mass):
        try:
            if len(p) == 3 and len(v) == 3 and type(mass) == float:
                self.pos = coords(*p)
                self.v = velocities(*v)
                self.m = mass
            else:
                raise Exception("Coordinates and velocities should be vectors of three floats each, mass needs to be float")
        except:
            raise Exception("position and velocity arguments should be vectors of length 3")

    def __repr__(self):
        posStr = self.pos.__repr__()
        vStr = self.v.__repr__()
        return "{0}\n{1}\nm = {2:.2f}".format(posStr,vStr,self.m)

class system:
    """
    class made up of all constituent gravitating bodies, contains nBodies with each body as described in body class.
    Optionally applies relavistic correction to gravity assuming one body (massiveInd) is much more massive than the others,
    with corrections calculated on smaller bodies that are "close enough" (closeInds). 
    Also keeps track of timestep system is evolved with and total evolved time of system.
    """
    def __init__(self,bodies,Δt=1.,massiveInd=0,closeInds=[]):
        try:
            self.nBodies = len(bodies)
            self.bodies = bodies
            self.massiveInd = massiveInd
            self.closeInds = closeInds
            self.Δt = Δt
            self.ΔtMod = 1
            self.T = 0. #total evolved time kept in units of yrs
        except:
            raise Exception("Bodies should be passed in as a list, massiveInd as an int, closeInds as a list, Δt as float")
    
    def getCoords(self): #get body all body coordinates in np arrays, easier to vectorize/operate on
        x = np.array([self.bodies[i].pos.x for i in range(self.nBodies)])
        y = np.array([self.bodies[i].pos.y for i in range(self.nBodies)])
        z = np.array([self.bodies[i].pos.z for i in range(self.nBodies)])

        vx = np.array([self.bodies[i].v.x for i in range(self.nBodies)])
        vy = np.array([self.bodies[i].v.y for i in range(self.nBodies)])
        vz = np.array([self.bodies[i].v.z for i in range(self.nBodies)])

        coords = np.array([x,y,z,vx,vy,vz])
        return coords
    
    def getMasses(self):
        masses = np.array([self.bodies[i].m for i in range(self.nBodies)])
        return masses

    def update(self,G=6.67408313131313e-11,yrErr=1e-10,getCoords=getCoords,getMasses=getMasses):
        coords = getCoords(self); masses = getMasses(self)
        Δt = self.Δt
        accept,rescale = acceptTimeStep(self,Δt/self.ΔtMod,yrErr,G) #try undoing timestep shrink
        
        if accept == True:
            self.ΔtMod = rescale
        else:
            print("Timestep could not be made small enough to satisfy error tolerance, continuing anyways")

        self.T += self.Δt/(365*24*3600)/rescale
        
        x,y,z,vx,vy,vz = nBodyStep(coords,masses,self.nBodies,self.Δt,self.massiveInd,self.closeInds,G) #from functions.py
        for n in range(self.nBodies):
            self.bodies[n] = body([x[n],y[n],z[n]],[vx[n],vy[n],vz[n]],self.bodies[n].m)

    def __repr__(self):
        return "System of nBodies = {0}\naccess each body's info with system.bodies\nCurrent Δt = {1:.2g} years, system has evolved {2:.2g} total years".format(self.nBodies,self.Δt/(365*24*3600)/self.ΔtMod,self.T/(365*24*3600))



def initTestBodies(nBodies=3,Δt=1.,relativity=True,closeInds=[1],massiveInd=0):
    """
    for testing purposes of classes, initializes a random system with three bodies
    """
    bodyList = []
    for n in range(nBodies):
        x,y,z = np.random.uniform(-100,100,3); vx,vy,vz = np.random.uniform(-10,10,3)
        bodyList.append(body([x,y,z],[vx,vy,vz],np.random.uniform(0,10)))

    return system(bodyList,Δt,massiveInd,closeInds) if relativity == True else system(bodyList,Δt,0,[])


