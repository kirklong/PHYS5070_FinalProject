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

def nBodyStep(bodies,nBodies,Δt,G=6.67408313131313e-11):
    x = np.array([bodies[i].pos.x for i in range(nBodies)])
    y = np.array([bodies[i].pos.y for i in range(nBodies)])
    z = np.array([bodies[i].pos.z for i in range(nBodies)])
    vx = np.array([bodies[i].v.x for i in range(nBodies)])
    vy = np.array([bodies[i].v.y for i in range(nBodies)])
    vz = np.array([bodies[i].v.z for i in range(nBodies)])
    masses = np.array([bodies[i].m for i in range(nBodies)])
    coords = np.array([x,y,z,vx,vy,vz])
    k1 = Δt*Δr(coords,masses,nBodies,G)
    k2 = Δt*Δr(coords+k1/2,masses,nBodies,G)
    k3 = Δt*Δr(coords+k2/2,masses,nBodies,G)
    k4 = Δt*Δr(coords+k3,masses,nBodies,G)
    coords += (k1 + 2.*k2 + 2.*k3 +k4)/6.
    return coords
