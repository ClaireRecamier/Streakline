STREAKLINE

This project uses the streakline method to trace the path of a stellar stream, given its globular cluster.

Based on Ana Bonaca's streakline code in C, this project uses the Chapel language instead for easier integration with parallel programming.

Streakline will be integrated into the ChplUltra project to create an orbital integrator that tracks a stellar stream in an evolving Ultralight Dark Matter potential.

METHOD OVERVIEW

Input:
galactic potential parameters (q1, q2, qz, halo radius)
globular cluster parameters (mass of gc, mass loss of gc, plummer radius, initial position, initial velocity)
simulation parameters (number of timesteps, number of stars ejected)
integrator details (leapfrog)

From the initial position, the globular cluster is first evolved backwards in time for the desired number of timesteps. From this position, the globular cluster is then evolved forward that number of timesteps, ejecting stars (every set number of timesteps) at the tidal radius of the globular cluster plus a random position offset and at the angular velocity of the globular cluster plus a random velocity offset. 


STRUCTURE:

1. P5
  P5.chpl
  a. P5-NFW
  b. P5-pointmass
  c. RadialOffsets


2. streakline-Chapel

3. streakline-original
