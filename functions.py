#!/usr/bin/env python
import numpy as np
import copy

def Δr(coords,masses,nBodies,massiveInd,closeInds,G,c=2.99792458e8):
    """
    Function we will numerically integrate, returns change in position and velocity after *Δt.
    Uses relativistic correction for bodies specified in closeInds assuming limit where m << M (M = masses[massiveInd]).
    
    params: coords: the coordinates of bodies in system, generated with system.getCoords();
            masses: the masses of the bodies in system, generated with system.getMasses();
            nBodies: the number of bodies in the system, accessed with system.nBodies;
            massiveInd: the index of the most massive body, to apply relativity to, accessed with system.massiveInd;
            closeInds: the indices of the smaller bodies to apply relativity to, accessed with system.closeInds;
            G: Gravitational constant
            c: speed of light (defaults to SI value)

    returns: Δ: the change in velocity and position, in the shape of coords
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
            sep = np.sqrt((x[n]-x[massiveInd])**2+(y[n]-y[massiveInd])**2+(z[n]-z[massiveInd])**2)
            v2 = (vx[n]-vx[massiveInd])**2 + (vy[n]-vy[massiveInd])**2 + (vz[n]-vz[massiveInd])**2
            # term1 = G*masses[massiveInd]/(c**2 * sep**3)
            # term2 = 4*G*masses[massiveInd]/sep
            # term3 = -v2
            # term4 = 4*((vx[n]-vx[massiveInd])*(x[n]-x[massiveInd])+(vy[n]-vy[massiveInd])*(y[n]-y[massiveInd])+(vz[n]-vz[massiveInd])*(z[n]-z[massiveInd]))
            # Δvx += term1*(term2*(x[n]-x[massiveInd])+term3*(x[n]-x[massiveInd])+term4*(vx[n]-vx[massiveInd]))
            # Δvy += term1*(term2*(y[n]-y[massiveInd])+term3*(y[n]-y[massiveInd])+term4*(vy[n]-vy[massiveInd]))
            # Δvz += term1*(term2*(z[n]-z[massiveInd])+term3*(z[n]-z[massiveInd])+term4*(vz[n]-vz[massiveInd]))
            
            
            Δvx = -G*masses[massiveInd]/(c**2*sep**3)*((x[n]-x[massiveInd])*(c**2 - 4*G*masses[massiveInd]/sep + v2) - 4*(vx[n]-vx[massiveInd])*((vx[n]-vx[massiveInd])*(x[n]-x[massiveInd])))
            Δvy = -G*masses[massiveInd]/(c**2*sep**3)*((y[n]-y[massiveInd])*(c**2 - 4*G*masses[massiveInd]/sep + v2) - 4*(vy[n]-vy[massiveInd])*((vy[n]-vy[massiveInd])*(y[n]-y[massiveInd])))
            Δvz = -G*masses[massiveInd]/(c**2*sep**3)*((z[n]-z[massiveInd])*(c**2 - 4*G*masses[massiveInd]/sep + v2) - 4*(vz[n]-vz[massiveInd])*((vz[n]-vz[massiveInd])*(z[n]-z[massiveInd])))


        Δ[3][n] = Δvx; Δ[4][n] = Δvy; Δ[5][n] = Δvz
        

    Δ[0] = vx; Δ[1] = vy; Δ[2] = vz
    return Δ

def nBodyStep(coords,masses,nBodies,Δt,massiveInd,closeInds,G=6.67408313131313e-11,integrator="RK4"): #integrating fx
    """
    RK4 and VV integrators using Δr

    params: coords: the coordinates of bodies in system, generated with system.getCoords();
            masses: the masses of the bodies in system, generated with system.getMasses();
            nBodies: the number of bodies in the system, accessed with system.nBodies;
            massiveInd: the index of the most massive body, to apply relativity to, accessed with system.massiveInd;
            closeInds: the indices of the smaller bodies to apply relativity to, accessed with system.closeInds;
            G: Gravitational constant (defaults to SI value)
            integrator: either "RK4" (4th order Runge-Kutta) or "VV" (Velocity-Verlet) integration scheme, defaults to RK4
    
    returns: new coords after update, of shape coords
    """
    if integrator == "RK4": #Runge-Kutta 4th order, explicit energy loss
        k1 = Δt*Δr(coords,masses,nBodies,massiveInd,closeInds,G)
        k2 = Δt*Δr(coords+k1/2,masses,nBodies,massiveInd,closeInds,G)
        k3 = Δt*Δr(coords+k2/2,masses,nBodies,massiveInd,closeInds,G)
        k4 = Δt*Δr(coords+k3,masses,nBodies,massiveInd,closeInds,G)
        coords += (k1 + 2.*k2 + 2.*k3 +k4)/6.
        return coords

    elif integrator == "VV": #Velocity-Verlet algorithm, symplectic but only second order
        Δ = Δr(coords,masses,nBodies,massiveInd,closeInds,G)
        coords[0:3] += Δ[0:3]*Δt #drift positions, v*dt
        coords[3:] += Δ[3:]*Δt/2 #kick velocities halfway
        Δ = Δr(coords,masses,nBodies,massiveInd,closeInds,G) #calculate again
        coords[3:] += Δ[3:]*Δt/2 #finish updating velocities from accelerations
        return coords

    else:
        raise Exception("Integrator not recognized: choose from VV (velocity-verlet) or RK4 (4th order Runge-Kutta)")

def E(coords,masses,nBodies,G=6.67408313131313e-11):
    """
    Calculates the total energy (gravitational potential + kinetic) of the system, does not include relativity
    
    params: coords: the coordinates of bodies in system, generated with system.getCoords();
            masses: the masses of the bodies in system, generated with system.getMasses();
            nBodies: the number of bodies in the system, accessed with system.nBodies;
            G: Gravitational constant, defaults to SI value
    
    returns: total energy of system
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

    params: system: object generated with system class; Δt: proposed timestep, i.e. system.Δt;
            yrErr: maximum acceptable energy loss per year, defaults to 1e-10; 
            G: Gravitational constant, defaults to SI value

    returns: accept: Boolean, whether we accept the timestep or not; rescale: how much we had to modify the timestep to keep within acceptable error
    """
    yearSec = 365*24*3600.
    maxErr = Δt/yearSec*yrErr
    Ei = E(system.getCoords(),system.getMasses(),system.nBodies,G)
    Ef = 0.; rescale = 1
    accept = False
    while accept == False and rescale <2**18:
        tmpCoords = copy.deepcopy(system.getCoords())
        coords = nBodyStep(tmpCoords,system.getMasses(),system.nBodies,Δt/rescale,system.massiveInd,system.closeInds,G)
        Ef = E(coords,system.getMasses(),system.nBodies,G)
        if np.abs((Ef-Ei)/Ei) < maxErr:
            accept = True
        else:
            rescale *= 2
    return accept, rescale





