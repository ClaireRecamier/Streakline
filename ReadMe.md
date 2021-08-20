STREAKLINE

This project uses the streakline method to trace the path of a stellar stream, given its progenitor globular cluster.

Based on Ana Bonaca's streakline code in C, this streakline uses the Chapel language instead for easier integration with parallel programming.

Streakline will work in tandem with the ChplUltra project to track a stellar stream in an evolving Ultralight Dark Matter potential.

METHOD OVERVIEW

Input:
a. galactic potential parameters (q1, q2, qz, halo radius)
b. globular cluster parameters (mass of gc, mass loss of gc, plummer radius, initial position, initial velocity)
c. simulation parameters (number of timesteps, number of stars ejected, integrator type)

From the initial position, the globular cluster is first evolved backwards in time for the desired number of timesteps. From this position, the globular cluster then moves forward that number of timesteps, ejecting pairs of stars every set number of timesteps. For each timestep, the position of the globular cluster is updated, followed by the position of all previously released stream particles. The program stores the position of the globular cluster at every timestep, but only stores the position of stream particles at the current timestep.

Pairs of stars are ejected every Mth timestep at the two lagrange points of the globular cluster, plus a random position offset, creating a leading tail of stars closer to the central potential and a trailing tail of stars further away from the central potential. Stars are ejected with a velocity equal to the angular velocity of the globular cluster plus a random velocity offset.

The orbit of the globular cluster is determined by the force due to the central potential, which can be either a pointmass or have an NFW density profile. Ejected stars, however, are subject to the force due to the central potential, and the force due to the potential of the globular cluster. In calculating the potential of the globular cluster, it is treated as a Plummer sphere.

STRUCTURE:

1. P5
  P5.chpl
  P5.c
  a. P5-NFW
      P5 NFW mathematica notebook
  b. P5-pointmass
      P5 pointmass mathematica notebook
  c. RadialOffsets
      Radial Offsets mathematica notebook


2. streakline-Chapel
  streakline.chpl

3. streakline-original
  streakline.c
