#!/usr/bin/env python
import numpy as np
def Δr(coords,masses,nBodies,G):
    x,y,z,vx,vy,vz = coords
    Δ = np.copy(coords)
    for n in range(nBodies):
        Δvx = 0.; Δvy = 0.; Δvz = 0.
        mask = np.arange(nBodies) != n
        sep = np.sqrt((x[n]-x[mask])**2+(y[n]-y[mask])**2+(z[n]-z[mask])**2)
        Δvx = np.sum(-G*masses[mask]*(x[n]-x[mask])/sep**3)
        Δvy = np.sum(-G*masses[mask]*(y[n]-y[mask])/sep**3)
        Δvz = np.sum(-G*masses[mask]*(z[n]-z[mask])/sep**3)
        Δ[3][n] = Δvx; Δ[4][n] = Δvy; Δ[5][n] = Δvz
    Δ[0] = vx; Δ[1] = vy; Δ[2] = vz
    return Δ

def nBodyStep(coords,masses,nBodies,Δt,G=6.67408313131313e-11): #RK4
    k1 = Δt*Δr(coords,masses,nBodies,G)
    k2 = Δt*Δr(coords+k1/2,masses,nBodies,G)
    k3 = Δt*Δr(coords+k2/2,masses,nBodies,G)
    k4 = Δt*Δr(coords+k3,masses,nBodies,G)
    coords += (k1 + 2.*k2 + 2.*k3 +k4)/6.
    return coords

def E(coords,masses,nBodies,G=6.67408313131313e-11):
    d = lambda coords,i,n : np.sqrt((coords[0][n]-coords[0][i])**2+(coords[1][n]-coords[1][i])**2+(coords[2][n]-coords[2][i])**2) 

    totalU = 0.
    x,y,z,vx,vy,vz = coords
    for n in range(nBodies-1):
        totalU += -G*np.sum(masses[n]*masses[n+1:]/d(coords,n,np.arange(nBodies-(n+1))+n+1))
    totalK = np.sum(masses*(vx**2+vy**2+vz**2))/2.
    return totalK + totalU


