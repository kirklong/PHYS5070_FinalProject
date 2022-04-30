#!/usr/bin/env python
from mimetypes import init
from classes import *
import pickle,sys
import matplotlib.pyplot as plt

def trackCompletion(place,stop,strLen): #percent output tracker
    string="Currently at T = {0:.2f} / {1:.2f} years".format(place,stop)
    sys.stdout.write("\r") #this "moves the cursor" to the beginning of the I0 line
    sys.stdout.write(" "*strLen) #this "clears" whatever was on the line last time by writing whitespace
    sys.stdout.write("\r") #move the cursor back to the start again
    sys.stdout.write(string) #display the current percent we are at
    sys.stdout.flush() #flush finishes call to print() (this is like what's under the hood of print function)
    strLen=len(string) #return the new string length for next function call
    return strLen

def getPlotData(system,cadence=1e-3,stopT=1,yrErr=1e-3,integrator="RK4"):
    strLen = 0; oldT = 0; place = 0
    x = np.zeros((system.nBodies,int(np.round(stopT/cadence))))
    y = np.zeros((system.nBodies,int(np.round(stopT/cadence))))
    z = np.zeros((system.nBodies,int(np.round(stopT/cadence))))
    t = np.zeros(int(np.round(stopT/cadence)))
    while system.T < stopT:
        system.update(yrErr=yrErr,integrator=integrator)
        #system.bodies[0].pos = coords(0.,0.,0.) #keep the Sun from moving
        #system.bodies[0].v = velocities(0.,0.,0.)
        if system.T-oldT > cadence:
            strLen = trackCompletion(system.T,stopT,strLen)
            oldT = system.T
            xtmp,ytmp,ztmp,vx,vy,vz = system.getCoords()
            x[:,place] = xtmp
            y[:,place] = ytmp
            z[:,place] = ztmp
            t[place] = system.T
            place += 1

    return system,t[0:place],x[:,0:place],y[:,0:place],z[:,0:place]

def toMFrame(xM,yM,zM):
    φ = 48.331/180*np.pi #deg to rad, longitude of ascending node
    θ = 7.005/180*np.pi #deg to rad, inclination relative to ecliptic
    ψ = 29.124/180*np.pi #deg to rad, argument of periapsis

    a11 = np.cos(ψ)*np.cos(φ)-np.cos(θ)*np.sin(φ)*np.sin(ψ)
    a12 = np.cos(ψ)*np.sin(φ)+np.cos(θ)*np.cos(φ)*np.sin(ψ)
    a13 = np.sin(ψ)*np.sin(θ)
    a21 = -np.sin(ψ)*np.cos(φ)-np.cos(θ)*np.sin(φ)*np.cos(ψ)
    a22 = -np.sin(ψ)*np.sin(φ)+np.cos(θ)*np.cos(φ)*np.cos(ψ)
    a23 = np.cos(ψ)*np.sin(θ)
    a31 = np.sin(θ)*np.sin(φ)
    a32 = -np.sin(θ)*np.cos(φ)
    a33 = np.cos(θ)

    #rotate orbit to be in orbital plane of Mercury instead of ecliptic
    xnew = a11*xM+a12*yM+a13*zM
    ynew = a21*xM+a22*yM+a23*zM
    znew = a31*xM+a32*yM+a33*zM #this should now be basically zero
    return xnew,ynew

Msun = 1.9885e30; AU=1.495978707e11
Mercury = body(np.array([1.467212043857380E+07,4.358719474022335E+07,2.216070698313603E+06])*1e3,np.array([-5.592391315445695E+01,1.739218372449630E+01,6.550993329084459E+00])*1e3,3.302e23)
Sun = body([0.,0.,0.,],[0.,0.,0.,],Msun) #sun-centered version

yearSec = 365*24*3600
MP = 87.96926/365*yearSec
SunMercury = system([Sun,Mercury],Δt=MP/100,massiveInd=0,closeInds=[1])
initialOrbit = getPlotData(SunMercury,stopT=MP/yearSec+MP/100/yearSec,cadence=MP/100/yearSec)
stop = 100*3750 #15,000 centuries = 180 degrees, so 3750 centuries ~ 45 degrees, should take ~ 24 hours (slow python code woohoo)
strlen = 0

while SunMercury.T < stop:
    strlen = trackCompletion(SunMercury.T,stop,strlen)
    SunMercury.update()
    
endOrbit = getPlotData(SunMercury,stopT=stop+MP/yearSec+MP/100/yearSec,cadence=MP/100/yearSec)

with open('orbits.pkl','wb') as f:
    pickle.dump([initialOrbit,endOrbit],f)

font = {'family' : 'DejaVu Serif',
    'weight' : 'normal',
    'size'   : 16}
plt.rc('font', **font) #set all plot attribute defaults

def pltFormatter(fig,axList,**kwargs):
    for ax in axList:
        ax.minorticks_on()
        ax.grid(b=True,which="major",alpha=0.5)
        ax.grid(b=True,which="minor",alpha=0.3)
        legend=kwargs.get("legend")
        if legend != None:
            legend.get_frame().set_edgecolor('black') 

fig = plt.figure(figsize=(10,10),facecolor='white')
ax = fig.gca()

xnew,ynew = toMFrame(*initialOrbit[2:])
ax.plot(xnew[0,:]/AU,ynew[0,:]/AU,label="Sun",marker=".",ms=20,c="gold")
ax.plot(xnew[1,:]/AU,ynew[1,:]/AU,label="starting orbit",marker=".",ms=4,c="forestgreen",alpha=0.5)
xnew,ynew = toMFrame(*endOrbit[2:])
ax.plot(xnew[1,:]/AU,ynew[1,:]/AU,label="final orbit",marker=".",ms=4,c="crimson",alpha=0.5)
ax.set_xlabel("X [AU]"); ax.set_ylabel("Y [AU]")
l = ax.legend()
pltFormatter(fig,[ax],legend=l)
ax.set_facecolor('k')
fig.tight_layout()
fig.savefig("orbits.png")
