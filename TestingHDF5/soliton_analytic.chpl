use IO;
use Math;
use Random;
//DIMENSIONLESS VERSION
//math constants
const pi = 3.141592653589793;
const Msun = 4 * pi**2; //AU
const mau = 6.68458e-12; //meters to AU
const kmkpc = 3.24078e-17;
const kpcau = 2.063 * (10**8); //kpc to au
const aukpc = 4.84814e-9;
const yrsec: real = 60*60*24*365.24; //years to seconds
const secyr: real = 1 / yrsec;
const kpckm = 1.0 / kmkpc;
const G = 6.67430e-11; //m^3 / (kg sec^2)
const G_code = 6.67430e-11 * ((1e-3)**3) * (yrsec)**2 * kmkpc**3 * (1/toCodeTime(1))**2  * toCodeLength(1)**3 * (1/toCodeMass(1)); //G in km^3 / (kg * sec^2)
writeln("G_code ",G_code);
const hbar = 1.0545718e-34; //((m^2)*(kg/s))
const parsec = 3.08571e16 ; // m
const H0 =(67.7e3)/(parsec*(1e6)) ; // 1/sec
const omegaM0 = 0.31 ;
const maxion = 1e-22 * 1.783 * 1e-36;
//setup constants
config const r0 = toCodeLength(15.0); // radius of ring, in km
config const integrator = 0;   //set Integrator to LF
config const N = 100000;//6000 set number of timesteps to equal total of 6 Gyr
var dt = toCodeTime(1000000.0); //set timestep to seconds per Myr
config const offsetOption = 0; //0 = no radial offsets, 1 = using C generated random numbers, 2 = generating own random numbers
config const Ne = 1; //number of particles in ring
config const nbins = 4;
var pot = 5; //set galactic potential to that of soliton
var period: real = 0.0;
config const po = (0.5**4) * 27492.260803351877;
config const rc = 2.0*0.05359134269304645;

//calculate potential at every point, output to file
//box ranges from x(-1.5,1.5) y(-1.5,1.5) z(-1.5,1.5)
proc load_pot(calcpar,pot) { //creates acceleration box in km/seconds squared
  var par: [0..5] real = [417.0, 36.54 * kpckm, 90.0 * pi / 180, 1.0, 1.0 ,1.0];
  if pot == 3 { //if using triaxial NFW potential
    //assuming par = [V, rhalo, phi, q_1, q_2, q_z]
    //calcpar = [GM, c1, c2, c3, c4, rhalo]
    var cosphi: real = cos(par[2]);
    var sinphi: real = sin(par[2]);
    calcpar[0] = par[0]*par[0]*par[1]; //GM in km cubed / seconds squared
    calcpar[1] = (cosphi**2)/(par[3]*par[3]) + (sinphi**2)/(par[4]*par[4]);
    calcpar[2] = (cosphi**2)/(par[4]*par[4]) + (sinphi**2)/(par[3]*par[3]);
    calcpar[3] = 2*sinphi*cosphi*(1/(par[3]**2) - 1/(par[4]**2));
    calcpar[4] = 1/(par[5]*par[5]);
    calcpar[5] = par[1];
    var acc: 3*real;
    //acc = force(r0 * pos[j],pot,calcpar); //convert pos to dimensions and get force
    writeln("GM in km cubed / seconds squared ", calcpar[0]);
    acc = force((r0,0.0,0.0),pot,calcpar);
    writeln("centripetal acceleration in km/sec^2 ",acc[0]," ",acc[1]," ",acc[2]);
    calcpar[5] = par[1] / r0; //Rhalo becomes dimensionless
    period = 2 * pi * r0 / sqrt(len(acc) * r0); //in seconds; 2pi * r/v

    //period = 2 * pi * r0 / magVel; //in seconds; 2pi * r/v
    calcpar[0] = calcpar[0] * (period ** 2)/ (r0 ** 3); //GM becomes dimensionless
    dt = dt / period; //divide dt by period to make it dimensionless
    writeln("period ", period);
    writeln("dt ", dt);

    }
  else if pot == 0 {
    calcpar[0] = G * 1000000000.0 * (1.989e30); //GM in km^3 / sec^2
  }
}

proc force(pos,pot,calcpar){
  var acc: 3*real;
  var r: real;
  if pot == 0 { //if using point mass potential
    r = len(pos);
    //acc = (-1,-1,-1)*(Msun * pos)/dist**3; //assumes the sun stays at origin
    //acc = (-1,-1,-1)*(100000000000 * Msun * pos)/r**3;
    //writeln("acc ",acc);
    acc = (-1,-1,-1)*(calcpar[0] * pos)/r**3;

  }
  else if pot == 3 { //NFW triaxial potential, C version
    //calcpar = [GM, c1, c2, c3, c4, rhalo]
    r = sqrt(calcpar[1]*pos[0]*pos[0] + calcpar[2]*pos[1]*pos[1] + calcpar[3]*pos[0]*pos[1] + calcpar[4]*pos[2]*pos[2]);
    var aux: real = 0.5 * calcpar[0] / (r**3) * (1.0/(1.0 + calcpar[5]/r)-log(1.0+r/calcpar[5]));

    acc[0]=aux*(2*calcpar[1]*pos[0] + calcpar[3]*pos[1]);
    acc[1]=aux*(2*calcpar[2]*pos[1] + calcpar[3]*pos[0]);
    acc[2]=aux*(2*calcpar[4]*pos[2]);
    //writeln("acc magnitude ", len(acc));
  }
  else if pot == 5 {
    //var pos1: 3*real = (pos[0]+0.001,pos[1]+0.001,pos[2]+0.001);
    r = len(pos);
    var height: real = 0.0;
    var dr: real = 0.1;
    var mass: real = 0.0;
    var bound: int = (r/dr): int;
    var a: real = 0.0;
    for i in 0..bound {
      const r1 = i * dr;
      height = po * 4 * pi * (r1**2) / (1 + 0.091 * (i * dr/rc)**2)**8; //
      mass = mass + (height * dr);
      //a += -G * height * dr / (i*dr)**2;
    }
    acc = (- G_code * mass / r**2) * (pos/r);
    //acc = a * (pos/r);
    //writeln("A ",(- G_code * mass / r**2));
  }
  return acc;
}

//initialize particles in a ring
proc init_ring (pos, vel, AM, SD, PS, calcpar,WritingChannel,AMWritingChannel,PSWritingChannel) {
  //assign counterclockwise from rightmost point on x axis.
  var posAngle, velAngle, magVel: real;
  var acc: 3*real;
  //acc = force((r0,0.0,0.0),pot,calcpar);
  acc = force((r0,0.0,0.0),pot,calcpar);
  //magVel = sqrt(len(acc) * r0) * period / r0;
  //magVel = sqrt(len(acc) * r0); //already dimensionless
  magVel = sqrt(len(acc)*r0);
  //writeln("magVel ",magVel);
  var period: real = 2 * pi * r0 / magVel; //in seconds; 2pi * r/v
  //writeln("dimensionless dt ", dt / period);

  //var velAngle: real;
  for j in 1..Ne { //Ne=number of desired particles
    posAngle = j * 2.0 * pi / Ne; //cycle thru angles in a circle

    pos[j][0] = r0 * cos(posAngle);
    pos[j][1] = r0 * sin(posAngle);
    pos[j][2] = 0.0;
    //pos[j][0] = cos(posAngle); //dimensionless
    //pos[j][1] = sin(posAngle);
    //pos[j][2] = 0.0;

    velAngle = (pi / 2.0) - posAngle; //angles are complementary

    vel[j][0] = magVel * cos(velAngle); //dimensionless
    vel[j][1] = -1.0 * magVel * sin(velAngle);
    vel[j][2] = 0.0;

    AM[j] = angMom(pos[j],vel[j]);
    //writeln("initial pos ", pos[j]);
    //writeln("initial vel ", vel[j]);
    WritingChannel.write(pos[j][0]*vel[j][0] + pos[j][1]*vel[j][1] + pos[j][2]*vel[j][2],",0.0,0.0,",len(vel[j]),",",len(pos[j]),",");
    WritingChannel.write(pos[j][0],",",pos[j][1],",",pos[j][2],",");
    WritingChannel.write(vel[j][0],",",vel[j][1],",",vel[j][2],"\n");
    AMWritingChannel.write(AM[j],",");
  }
  pow_spec(pos,vel,PS,PSWritingChannel);
  AMWritingChannel.write("\n");

  //iterate through 2pi radians divided by 10 (Ne) particles.
  //for each, multiply r0 by sin to get y, r0 by cosine to get x
  //calculate magnitude of centripetal velocity vector for each particles
  //angle velocity vector makes with axis is 90 - theta.
  //multiply magnitude of velocity vector with sin and cosine of angle.
}

//orbit procedure: advances cluster in position and velocity using integrator of choice by N timesteps
proc fwd_orbit (pos, vel, AM, SD, PS, pot, integrator, N, dt, calcpar,WritingChannel, SDWritingChannel,AMWritingChannel,PSWritingChannel) {
    if integrator == 0  { //if leapfrog
      //move velocity of all particles a halfstep forward

      for j in 1..Ne {
        halfstep(pos[j],vel[j],pot,dt,1.0, calcpar);
        AM[j] = angMom(pos[j],vel[j]);
        //writeln("vel after halfstep ",vel[j]);
        //writeln("pos after halfstep ",pos[j]);

        //WritingChannel.write(pos[j][0] * aukpc,",",pos[j][1] * aukpc,",",pos[j][2] * aukpc,",");
        //WritingChannel.write(vel[j][0] * aukpc,",",vel[j][1] * aukpc,",",vel[j][2] * aukpc,"\n");
        //AMWritingChannel.write(AM[j],",");
      }
      writeln("first halfstep forward complete");
      //AMWritingChannel.write("\n");
      for i in 0..N { //for each subsequent timestep
        //move velocity of all particles a fullstep forwards
        for j in 1..Ne {
          stream_step(pos[j], vel[j], dt, calcpar,WritingChannel);
          //writeln("streamstep complete");
          AM[j] = angMom(pos[j],vel[j]);
          //writing position and velocity of to file (one particle per line)
          WritingChannel.write(pos[j][0],",",pos[j][1],",",pos[j][2],",");
          WritingChannel.write(vel[j][0],",",vel[j][1],",",vel[j][2],"\n");
          //writing angular momentum to a diff file
          AMWritingChannel.write(AM[j],",");
          }
        AMWritingChannel.write("\n");
        SD[i] = stdev(AM);
        pow_spec(pos,vel,PS,PSWritingChannel);
        //writeln("fullstep complete");
        /* //print snapshot at 5000th timestep
        if i == N-1000 {
          for j in 1..k {
            myWritingChannel.write(pos_lead[j][0],",",pos_lead[j][1],",",vel_lead[j][0],",",vel_lead[j][1],",",pos_trail[j][0],",",pos_trail[j][1],",",vel_trail[j][0],",",vel_trail[j][1],"\n");
          }
        }
        */
      }
      //for final timestep, move everything back
      for j in 1..Ne {
        halfstep(pos[j],vel[j],pot,dt,-1.0, calcpar);
      }
      for i in 0..N {
        SDWritingChannel.write(i,",",SD[i] * (r0/period),"\n");
      }
  }
}

//procedure to advance ejected particles by a timestep. only subject to force from central pot
proc stream_step(ref pos, ref vel, dt, calcpar,WritingChannel) {
  //writeln("pos ",pos);
  //update position of jth particle
  pos += dt * vel;
  //writeln("vel ",vel);

  //writeln("pos ",pos);

  //search acceleration
  var a: 3*real;
  //writeln("searching for force");
  a = force(pos,pot,calcpar); //acc made dimensionless
  //writeln("found force");

  //writeln("acceleration ",a);
  //update velocity of jth particle
  vel += dt * a;
  //writeln(vel);

  //writeln("position ",pos);
  //writeln("velocity ", vel);
  WritingChannel.write(pos[0]*vel[0] + pos[1]*vel[1] + pos[2]*vel[2],",");
  WritingChannel.write(a[0]*vel[0] + a[1]*vel[1] + a[2]*vel[2],",");
  WritingChannel.write(len(a),",");
  WritingChannel.write(len(vel),",");
  WritingChannel.write(len(pos),",");
  //writeln("magnitude velocity: ",len(vel));
  //writeln("magnitude position: ",len(pos));
  //writeln(a[0]*vel[0] + a[1]*vel[1] + a[2]*vel[2]);
}

//shifts velocity by a halfstep in direction of sign
proc halfstep(p, ref v, pot, dt, sign, calcpar) {
  //writeln("pos halfstep",p);
  //writeln("vel halfstep",v);
  var signed_dt: real = sign * dt;
  var a: 3*real;
  a = force(p,pot,calcpar);
  //writeln("acc halfstep ",a);
  //writeln("dot force and vel ",a[0]*v[0] + a[1]*v[1] + a[2]*v[2]);

  //writeln("signed dt ",signed_dt," calc from halfstep: ",0.5 * signed_dt * a);
  //writeln("a from halfstep ",a);
  v += 0.5 * signed_dt * a;
  //writeln("vel halfstep",v);


}

proc angMom(pos,vel){
  var omega: 3*real;
  omega = cross(pos, vel); //r x v
  //omega = omega / (len(pos)**2); //(r x v)/r^2
  return len(omega);
}
proc cross((x1,y1,z1),(x2,y2,z2)) {
  var cross_product: 3*real;
  cross_product[0] = (y1 * z2) - (z1 * y2);
  cross_product[1] = (z1 * x2) - (x1 * z2);
  cross_product[2] = (x1 * y2) - (y1 * x2);
  return cross_product;
}
proc stdev(AM){
  var mean, SD: real = 0.0;
  for j in 1..Ne { //iterate over all stream particles
    mean += AM[j];
  }
  mean = mean / Ne;
  for j in 1..Ne {
    SD += (AM[j] - mean)**2;
  }
  SD = SD / Ne;
  SD = sqrt(SD);
  //return SD/r0;
  return SD;
}
proc len((x,y,z)) {
  return sqrt(x**2 + y**2 + z**2);
}
proc pow_spec (pos,vel, PS,PSWritingChannel) {//gets computed at each timstep
  //keep track of index of particle
  var j: int = 1;
  var avg: real = Ne / nbins;
  //loop through number of bins, aka sections of angle
  var low_bound, up_bound, posAngle: real = 0.0;
  posAngle = acos( pos[j][0] / len(pos[j]) );
  for i in 1..nbins{
    low_bound = (i - 1.0) * (2.0 * pi / nbins);
    up_bound = i * (2.0 * pi / nbins);
    //writeln("low: ",low_bound);
    //writeln("high ", up_bound);
    //calculate angle of position of current particle
    //while angle of position is in current bin, iterate over particles
    PS[i] = 0;
    while posAngle > low_bound && posAngle <= up_bound{
      //writeln("j ", j);
      //writeln(pos[j]);
      //writeln(posAngle);
      PS[i] += 1; //add one to current bin
      j += 1; //move on to next particle
      posAngle = acos( pos[j][0] / len(pos[j]) ); //update angle
      if pos[j][1] < 0{  //if y is neg, aka in lower quad,
        posAngle = posAngle + 2 * (pi - posAngle);
        }
    }
    //PSWritingChannel.write(i,",",PS[i] / avg,",");
    PSWritingChannel.write(PS[i] / avg,",");

  }
  PSWritingChannel.write("\n");
}

proc toCodeTime(yrs){
  return yrs / (75.4212 * 10**9);
}
proc toCodeLength(kpc){
  return kpc / 38.3609;
}

proc toCodeMass(kg) {
  return kg/4.42837e36;
}

proc main(){
  //create array of tuples to hold position and velocity and velocity at each timestep

  var pos: [1..Ne] 3*real;
  var vel: [1..Ne] 3*real;
  var AM: [1..Ne] real; //stores angular momentum of each particle, updated at each timestep
  var SD: [0..N] real; //stores std dev at each timestep
  var PS: [1..nbins] real; //whole array gets updated at each timestep
  //hardcode galactic potential parameters
  var calcpar: [0..5] real;

  var chfile = open("../TestingHDF5/mathematica/bp15.csv",iomode.cw); //create test.csv and open
  var SDfile = open("../TestingHDF5/mathematica/bpSD15.csv",iomode.cw); //create test.csv and open
  var AMfile = open("../TestingHDF5/mathematica/bpAM15.csv",iomode.cw); //create test.csv and open
  var PSfile = open("../TestingHDF5/mathematica/bpPS15.csv",iomode.cw); //create test.csv and open
  var PSWritingChannel = PSfile.writer(); //open writing channel to test.csv
  var SDWritingChannel = SDfile.writer(); //open writing channel to test.csv
  var WritingChannel = chfile.writer(); //open writing channel to test.csv
  var AMWritingChannel = AMfile.writer(); //open writing channel to test.csv

  load_pot(calcpar,pot);
  writeln("loaded pot");

  //initialize initial position and velocity of every particle in ring
  init_ring(pos, vel, AM, SD, PS, calcpar, WritingChannel,AMWritingChannel,PSWritingChannel);
  //dt = dt / period;
  writeln("initialized ring");
  writeln("pos 1 ",pos[1]);
  writeln("vel 1", vel[1]);

  //
  /* writeln(1.0-0.015625*17/8); */
  /* writeln(force((-0.5,-0.5,-0.5),pot,calcpar)); */

  fwd_orbit(pos, vel, AM, SD, PS, pot, integrator, N, dt, calcpar,WritingChannel,SDWritingChannel,AMWritingChannel,PSWritingChannel);
  WritingChannel.close();
  SDWritingChannel.close();
  AMWritingChannel.close();
  PSWritingChannel.close();
  PSfile.close();
  AMfile.close();
  SDfile.close();
  chfile.close();
}
