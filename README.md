# PHYS 5070 Final Project &mdash; An n-body simulation of the solar system to verify the total precession of Mercury.

#### Author: Kirk Long


## The idea:

I've built a toy n-body simulator pretty much entirely from scratch that I'm using to simulate the solar system. My (perhaps overly optimistic) idea was that I'd be able to recover the precession of Mercury's orbit from the influence of the outer planets as well as general relativity. 

For reference, the observed total precession of Mercury's perihelion is ~5.75"/year. Newtonian/Keplerian dynamics from the outer planets account for 5.32"/year of this total value, but the remaining 0.43"/year was not explained until the advent of general relativity! I hope to recover both of these contributions. 

## The process: 

I decided to try a mostly object oriented approach for this project as I don't usually do that (usually I prefer to work in Julia, which is more functional). The file (`classes.py`)[classes.py] contains the functions that build these objects, culminating in the `system` class, which contains the entire "solar system" in my code. The `system` class has an `update` method that evolves the system forward in time according to a specified time step, with the option to choose between RK4 and a velocity-verlet integration scheme.

The file (`functions.py`)[functions.py] contains the physics of the simulation, most importantly in the `Δr` function. This is what I numerically integrate, with force interactions between Newtonian/Keplerian bodies going as usual like 1/r^2. Relativity can also be "turned on" for this function for specified bodies of the system (in my case I turn on relativity for the Mercury-Sun interaction only) using a weak-field approximation of GR, where the mass of the orbiting body (i.e. a planet like Mercury) is much less than the mass of the central object (i.e. a star like our Sun).

In this case I use the following approximation ((equation 2 from Parsa+2017)[https://www.eso.org/public/archives/announcements/pdf/ann17051a.pdf], but also given in many other places):

!(GR)[GREquation.png]

The first term is (as expected) the Newtonian form, whereas all the other terms are perturbations due to the curvature of space-time in GR caused by the Sun's large mass.

To simulate the solar system I first start with initial conditions from JPL in the ecliptic reference frame, queried via emailed jobs like: 

```
!$$SOF
EMAIL_ADDR=''
START_TIME = '2022-Apr-13 17:30:58'
STOP_TIME = '2022-Apr-13 17:30:59'
TABLE_TYPE = 'Vector'
REF_PLANE = 'Ecliptic'
CENTER = '@010'
COMMAND='301'
!$$EOF
```

Where `command` is the solar system body we want to retrieve data for (100 = Mercury, 200 = Venus, 300 = Earth, 301 = the Moon, etc.). This returns the most up-to-date ephemeride data (and much more), of which we use initial x,y,z,vx,vy,vz parameters. 

After inputting these initial parameters I then evolve the solar system in time. After the time evolution has finished I use Mercury's known longitude of ascending node, argument of perihelion, and inclination from the ecliptic to transform the coordinates into Mercury's orbital frame. I then find the point of periapsis for each orbit by calculating the set of minimum distances from the Sun in this frame, and I return the angle φ that corresponds to each. We can then plot this angle over time to see how it evolves, as this indicates precession! Since the precession is quite small over a single orbital period we have to integrate for long time scales to resolve this, and there is some inherent numerical noise related to the timestep as we will not always exactly sample where perihelion occurs. 

To calculate the total precession of Mercury we complete the above procedure using the "full" solar system, and to calculate just the effects from GR we simulate just the Sun-Mercury system.

## Tests:
Validation tests were performed on both integration schemes for the following cases:

1. Calculating the Earth's orbit in just a two-body interaction without relativity, and confirming that the periodicity is one year.
2. Calculating the Earth's orbit in the full n-body simulation without relativity, and again conforming that the periodicity is one year.
3. Calculating the Earth's orbit in just a two-body interaction with relativity turned on, and confirming that this does not affect the orbit.
4. Calculating Mercury's orbit in just a two-body interaction with relativity turned on, and confirming that this *does* change the orbit.

These tests + resulting plots can be found in the (`tests.ipynb`)[tests.ipynb] notebook.

## Results: 



## Conclusions: 

