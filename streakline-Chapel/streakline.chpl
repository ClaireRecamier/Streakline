const pi = 3.141592653589793;
const Msun = 4 * pi**2; //AU
config const integrator = 0;   //set Integrator to LF
config const N = 1;   //set number of timesteps
config const dt = 0.5; //set timestep
config const mcl = 0.5 * Msun;
config const M = 100.0; //particles are released every Mth timestep
var pot = 0; //set potential to that of pointmass
var Ne = ceil(N/M) : int; //number of particles released
var k: int = 0; //record how many particles have been released



proc main () {
  //create array of tuples to hold position and velocity and velocity at each timestep
  var pos: [0..N] 3*real;
  var vel: [0..N] 3*real;
  //hardcode initial position and velocity
  pos[0]=(10,0,0);
  vel[0]=(0,sqrt(Msun/len(pos[0])),0); //set velocity equal to centripetal velocity
  //vel[0] = (0,20,0);

  var pos_lead: [1..Ne] 3*real;
  var pos_trail: [1..Ne] 3*real;
  var vel_lead: [1..Ne] 3*real;
  var vel_trail: [1..Ne] 3*real;

  //hardcode values for stream particles
  k = 1;
  pos_lead[1] = (9,0,0);
  vel_lead[1] = (0,sqrt(Msun/len(pos_lead[1])),0);
  pos_trail[1] = (11,0,0);
  vel_trail [1] = (0,sqrt(Msun/len(pos_trail[1])),0);


  //writeln("initial energy ", energy(pos,vel,0));
  orbit(pos, vel, pot, integrator, N, dt, pos_lead, pos_trail, vel_lead, vel_trail);
  //writeln("final energy ", energy(pos,vel,N));
/*
  for i in 0..N do {
    //writeln(pos[i][0]); //print x coordinates
    //writeln(pos[i][1]); //print y coordinates
  } */

}

//orbit procedure: advances cluster in position and velocity using integrator of choice by N timesteps
proc orbit (pos, vel, pot, integrator, N, dt, pos_lead, pos_trail, vel_lead, vel_trail) {
  //should i store halfstep at index 0 or index 1
  if integrator == 0  { //if leapfrog
    //move velocity forward half a timestep
    halfstep(pos[0],vel[0],pot,dt,1.0);
    //writeln(pos_lead[1][0]);
    //make N full steps in pos and vel forwards
    //writeln(pos_lead[1][0]);
    for i in 1..N {//for each timestep
      //decrease mass

      leapfrog(pos,vel,i,dt);

      for j in 1..k {

        stream_step(pos_lead[j], vel_lead[j], pos[i], dt);
        //writeln(pos_lead[j][0]);
        stream_step(pos_trail[j], vel_trail[j], pos[i], dt);
      }

    }
    //move velocity backward half a timestep
    halfstep(pos[N],vel[N],pot,dt,-1.0);
  }
  /*
  else { //if RK
    for i in 1..N do
      RK();
  }
  */
}

//procedure to advance ejected particles by a timestep
proc stream_step(ref pos, ref vel, pos_cl, dt) {
  //update position of jth particle
  pos += dt * vel;
  //calculate acceleration from Msun
  var a: 3*real;
  a = force(pos,pot);
  //calculate plummer acceleration
  //a += force_plummer();
  //update velocity of jth particle
  vel += dt * a;
}

//procedure to eject particles
proc eject() {


}

//shifts velocity by a halfstep in direction of sign
proc halfstep(p,ref v,pot,dt,sign) {
  var signed_dt: real = sign * dt;
  var a: 3*real;
  a = force(p,pot);
  v += 0.5 * signed_dt * a;
}

//shifts velocity and position by a full step forward
proc leapfrog(pos, vel, i, dt) {
  var a: 3*real;
  pos[i] = pos[i-1] + dt * vel[i-1];
  a = force(pos[i],pot);
  vel[i] = vel[i-1] + dt * a;
  /*
  a = force(pos[i-1],pot);
  vel[i] = vel[i-1] + dt * a;
  pos[i] = pos[i-1] + dt * vel[i];
  */
}

proc RK(pos,vel,i,dt) {

}

proc force(pos,pot){
  var acc: 3*real;
  if pot == 0 { //if using point mass potential
    var dist: real = len(pos);
    acc = (-1,-1,-1)*(Msun * pos)/dist**3; //assumes the sun stays at origin
  }
  return acc;
}

proc force_plummer() {

}

proc len((x,y,z)) {
  return sqrt(x**2 + y**2 + z**2);
}

proc energy(pos,vel,i){
  var energy: real = 0.0;

  energy += 0.5 * mcl * len(vel[i])**2;
  energy -= Msun * mcl / len(pos[i]);

  return energy;
}
