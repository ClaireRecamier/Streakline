config const integrator = 0;   //set Integrator
config const N = 3;   //set number of timesteps
config const dt = 0.5; //set timestep

proc main () {
  //create array of tuples to hold position and velocity and velocity at each timestep
  var pos: [0..N] 3*real;
  var vel: [0..N] 3*real;
  //hardcode initial position and velocity
  pos[0]=[5,5,5];
  vel[0]=[1,1,1];
  var pot = 0; //set potential to that of pointmass
  orbit(pos, vel,pot, integrator,N,dt);

}

//orbit procedure takes in: initial pos and vel of cluster, output arrays, which potential, which integrator, number of timesteps,length of timestep, and direction
proc orbit (pos,vel, pot, integrator, N, dt) {
  //should i store halfstep at index 0 or index 1
  if integrator == 0  { //if leapfrog
    //move velocity forward half a timestep
    var v: [3] real = vel[0];
    var p: [3] real = pos[0];
    halfstep(p,v,pot,dt,1.0);
    vel[0] = v; //??
    //make N full steps in pos and vel forwards
    for i in 1..N do
      leapfrog(pos,vel,i,dt);
    //move velocity backward half a timestep
    v = vel[N];
    p = pos[N];
    halfstep(p,v,pot,dt,-1.0);
    vel[N] = v; //?
  }
  else if integrator == 1 { //if RK
    for i in 1..N do
      RK();
  }
}

//shifts velocity by a halfstep in direction of sign
proc halfstep(p,v,pot,dt,sign) {
  var signed_dt: real = sign * dt;
  var a: [3] real;
  force(p, a, pot);
  v += 0.5 * signed_dt * a;
}

//shifts velocity and position by a full step forward
proc leapfrog(pos, vel, i, dt) {
  var a: [3] real;
  var p: [3] real = pos[i];
  force(p,a,pot);
  vel[i] = vel[i-1] + dt * a;
  pos[i] = pos[i-1] + dt * vel[i];
}

proc RK(pos,vel,i,dt) {

}

proc force(pos,acc,pot){
  if pot == 0 { //if using point mass potential
    var dist: real = len(pos);
    acc = (G * Msun * pos)/dist**3;
  }

}

proc len((x,y,z)) {
  return sqrt(x**2 + y**2 + z**2);
}
