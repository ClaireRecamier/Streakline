use IO;

const pi = 3.141592653589793;
const Msun = 4 * pi**2; //AU
config const integrator = 0;   //set Integrator to LF
config const N = 1000;   //set number of timesteps
config const dt = 0.01; //set timestep
config const mcli: real = 9.54791938424326609e-04 * Msun; //initial mass of cluster
config const mclf: real = 9.53e-04 * Msun; //final mass of cluster
config const M = 200; //particles are released every Mth timestep
config const Rcl = 0.25; //radius of plummer core
var pot = 0; //set potential to that of pointmass
var Ne = ceil(N/M) : int; //number of particles released CHECK THIS
var k: int = 0; //record how many particles have been released
var dm = (mcli - mclf)/N; //amount of mass released per timesteps
var mcl: real = mcli; //current mass of cluster
var cfile = open("test3.csv",iomode.cw); //create test.csv and open
var myWritingChannel = cfile.writer(); //open writing channel to test.csv


proc main () {
  //writeln(force_plummer((1.0,0.0,0.0), 0.25));
  //create array of tuples to hold position and velocity and velocity at each timestep
  var pos: [0..N] 3*real;
  var vel: [0..N] 3*real;
  //hardcode initial position and velocity
  pos[0]=( 4.95008,0.0,0.0);
  vel[0]=( 0.0,2.89362,0.0); //set velocity equal to centripetal velocity

  var pos_lead: [1..Ne] 3*real;
  var pos_trail: [1..Ne] 3*real;
  var vel_lead: [1..Ne] 3*real;
  var vel_trail: [1..Ne] 3*real;

  /*
  //hardcode values for stream particles
  k = 1;
  pos_lead[1] = (9,0,0);
  vel_lead[1] = (0,sqrt(Msun/len(pos_lead[1])),0);
  pos_trail[1] = (11,0,0);
  vel_trail [1] = (0,sqrt(Msun/len(pos_trail[1])),0);
*/

  myWritingChannel.write("x cluster,y cluster,x cluster vel,y cluster vel,x lead trail,y lead tail,x vel lead tail, y vel lead tail,x trail tail,y trail tail,x vel trail tail, y vel trail tail\n");

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
    if integrator == 0  { //if leapfrog
      //move velocity forward half a timestep
      halfstep(pos[0],vel[0],pot,dt,1.0);
      myWritingChannel.write(pos[0][0],",",pos[0][1],",",vel[0][0],",",vel[0][1],"\n");
      //writeln("pos cluster after first halfstep: ",pos[0]);
      //writeln("vel cluster after first halfstep: ",vel[0]);
    for i in 1..N {//make N full steps in pos and vel forwards
      mcl -= dm;//decrease mass

      leapfrog(pos,vel,i,dt);
      myWritingChannel.write(pos[i][0],",",pos[i][1],",",vel[i][0],",",vel[i][1]);

      for j in 1..k {
        stream_step(pos_lead[j], vel_lead[j], pos[i], dt);
        stream_step(pos_trail[j], vel_trail[j], pos[i], dt);

        myWritingChannel.write(",",pos_lead[j][0],",",pos_lead[j][1],",",vel_lead[j][0],",",vel_lead[j][1],",",pos_trail[j][0],",",pos_trail[j][1],",",vel_trail[j][0],",",vel_trail[j][1]);
      }

      if i % M == 0 {
        k+=1;
        eject(pos[i],vel[i],pos_lead[k], vel_lead[k], pos_trail[k], vel_trail[k]);
        myWritingChannel.write(",",pos_lead[k][0],",",pos_lead[k][1],",",vel_lead[k][0],",",vel_lead[k][1],",",pos_trail[k][0],",",pos_trail[k][1],",",vel_trail[k][0],",",vel_trail[k][1]);
        //writeln("after ejecting ", vel_lead[k]);
        //writeln(pos_trail[k]);

      }

      myWritingChannel.write("\n");

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
  a += force_plummer(pos_cl - pos, Rcl);
  //update velocity of jth particle
  vel += dt * a;
}

//procedure to eject particles
proc eject(pos_cl, vel_cl, ref pos_lead, ref vel_lead, ref pos_trail, ref vel_trail) {
  //calculate angular velocity of Cluster
  var omega: 3*real;
  omega = cross(pos_cl, vel_cl); //r x v
  var r = len(pos_cl);
  omega = omega / (r**2); //(r x v)/r^2
  var om: real = len(omega);

  var Rj: 3*real = tidal_radius(pos_cl, vel_cl, om);//calculate tidal radius
  //initial position of particle is position of cluster plus or minus tidal radius
  writeln("Rj: ",Rj);
  //pos_lead = pos_cl - (Rj * (pos_cl/r)); //multiplied by unit vector
  //pos_trail = pos_cl + (Rj * (pos_cl/r));
  pos_lead = pos_cl - Rj;
  pos_trail = pos_cl + Rj;

  //initial velocity of particle is angular velocity of cluster times pos cluster - tidal radius
  //vel_lead = cross(omega, pos_lead); // v = omega cross r
  //vel_trail = cross(omega, pos_trail);

  var mag_vcl: real = len(vel_cl);
  vel_lead = vel_cl * (mag_vcl - len(Rj)* om)/(mag_vcl);
  vel_trail = vel_cl * (mag_vcl + len(Rj)* om)/mag_vcl;
  //writeln("cluster position at ejection ",pos_cl);
  //writeln("cluster velocity at ejection ",vel_cl);
  //writeln("pos_lead at ejection ",pos_lead);
  //writeln("vel_lead at ejection", vel_lead);
  //writeln("pos_trail at ejection ",pos_trail);
  //writeln("vel_trail at ejection", vel_trail);
}

proc cross((x1,y1,z1),(x2,y2,z2)) {
  var cross_product: 3*real;
  cross_product[0] = (y1 * z2) - (z1 * y2);
  cross_product[1] = (z1 * x2) - (x1 * z2);
  cross_product[2] = (x1 * y2) - (y1 * x2);
  return cross_product;
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

proc force_plummer(r, Rcl) {
  // dist from cluster to particle, cluster radius
  var dist = len(r);
  var raux: real = sqrt(dist**2 + Rcl**2);
  var acc: 3*real;
  acc = (mcl / raux**3 ) * r;
  return acc;
}

proc tidal_radius(pos_cl,vel_cl,omega) {
  //writeln("pos cl ",pos_cl," vel cl ",vel_cl," omega ",omega);
  var kpc: real = 2.063e+8; //kpc in AU
  var un_vec: 3*real = pos_cl/len(pos_cl);
  var delta = 0.02 * kpc * un_vec;
  var x1: 3*real = pos_cl - delta;
  //writeln("x1 " , x1);
  var x2: 3*real = pos_cl + delta;
  //writeln("x2 " , x2);
  var a1: 3*real = force(x1,pot);
  var a2: 3*real = force(x2,pot);
  //writeln("a1 " , a1);
  //writeln("a2 " , a2);

  var dpot: real = (len(a1) - len(a2))/(len(x1-x2));
  //var dpot: 3*real = (a1 - a2)/(len(x1-x2));
  //writeln("len x1-x2 ",len(x1-x2));
  //writeln("a1-a2 ",a1-a2);
  writeln("mcl ",mcl);
  return un_vec*((mcl/abs(omega*omega + dpot))**(1.0/3.0));
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
