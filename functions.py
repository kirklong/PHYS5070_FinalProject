#!/usr/bin/env python
import numpy as np
import copy

def Δr(coords,masses,nBodies,massiveInd,closeInds,G,c=2.99792458e8):
    """
    Function we will numerically integrate, returns change in position and velocity after *Δt.
    Uses relativistic correction for bodies specified in closeInds assuming limit where m << M (M = masses[massiveInd]).
    """
    x,y,z,vx,vy,vz = coords
    Δ = np.copy(coords)
    for n in range(nBodies):
        Δvx = 0.; Δvy = 0.; Δvz = 0.
        mask = np.arange(nBodies) != n
        sep = np.sqrt((x[n]-x[mask])**2+(y[n]-y[mask])**2+(z[n]-z[mask])**2)
        Δvx = np.sum(-G*masses[mask]*(x[n]-x[mask])/sep**3)
        Δvy = np.sum(-G*masses[mask]*(y[n]-y[mask])/sep**3)
        Δvz = np.sum(-G*masses[mask]*(z[n]-z[mask])/sep**3)
        if n in closeInds: #do relativistic correction assuming m << M, m = masses[closeInd] M = masses[massiveInd]
            sep = np.sqrt((x[massiveInd]-x[n])**2+(y[massiveInd]-y[n])**2+(z[massiveInd]-z[n])**2)
            α = G*masses[massiveInd]/sep**3/(c)**2
            β = 4*G*masses[n]/sep - ((vx[massiveInd]-vx[n])**2 + (vy[massiveInd]-vy[n])**2 + (vz[massiveInd]-vz[n]**2)) 
            γ = 4*((x[n]-x[massiveInd])*(vx[n]-vx[massiveInd]) + (y[n]-y[massiveInd])*(vy[n]-vy[massiveInd]) + (z[n]-z[massiveInd])*(vz[n]-vz[massiveInd]))
            Δvx += α*(β*(x[n]-x[massiveInd] + γ*(vx[n]-vx[massiveInd])))
            Δvy += α*(β*(y[n]-y[massiveInd] + γ*(vy[n]-vy[massiveInd])))
            Δvz += α*(β*(z[n]-z[massiveInd] + γ*(vz[n]-vz[massiveInd])))

        Δ[3][n] = Δvx; Δ[4][n] = Δvy; Δ[5][n] = Δvz
        

    Δ[0] = vx; Δ[1] = vy; Δ[2] = vz
    return Δ

def nBodyStep(coords,masses,nBodies,Δt,massiveInd,closeInds,G=6.67408313131313e-11): #RK4
    """
    RK4 integrator using Δr
    """
    k1 = Δt*Δr(coords,masses,nBodies,massiveInd,closeInds,G)
    k2 = Δt*Δr(coords+k1/2,masses,nBodies,massiveInd,closeInds,G)
    k3 = Δt*Δr(coords+k2/2,masses,nBodies,massiveInd,closeInds,G)
    k4 = Δt*Δr(coords+k3,masses,nBodies,massiveInd,closeInds,G)
    coords += (k1 + 2.*k2 + 2.*k3 +k4)/6.
    return coords

def E(coords,masses,nBodies,G=6.67408313131313e-11):
    """
    Calculates the total energy (gravitational potential + kinetic) of the system, does not include relativity
    """

    d = lambda coords,i,n : np.sqrt((coords[0][n]-coords[0][i])**2+(coords[1][n]-coords[1][i])**2+(coords[2][n]-coords[2][i])**2) 

    totalU = 0.
    x,y,z,vx,vy,vz = coords
    for n in range(nBodies-1):
        totalU += -G*np.sum(masses[n]*masses[n+1:]/d(coords,n,np.arange(nBodies-(n+1))+n+1))
    totalK = np.sum(masses*(vx**2+vy**2+vz**2))/2.
    return totalK + totalU

def acceptTimeStep(system,Δt,yrErr=1e-10,G=6.67408313131313e-11): #adaptive timestep
    """
    Determines if proposed timestep is accurate enough given acceptable numerical error,
    and if not attempts to shorten timestep until error is beneath threshold
    """
    yearSec = 365*24*3600.
    maxErr = Δt/yearSec*yrErr
    Ei = E(system.getCoords(),system.getMasses(),system.nBodies,G)
    Ef = 0.; rescale = 1
    accept = False
    while accept == False and rescale <2**10:
        tmpCoords = copy.deepcopy(system.getCoords())
        coords = nBodyStep(tmpCoords,system.getMasses(),system.nBodies,Δt/rescale,system.massiveInd,system.closeInds,G)
        Ef = E(coords,system.getMasses(),system.nBodies,G)
        if np.abs((Ef-Ei)/Ei) < maxErr:
            accept = True
        else:
            rescale *= 2
    return accept, rescale






