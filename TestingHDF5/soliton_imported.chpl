use IO;
use Math;
use Random;
use HDF5;
use H5Helper;
//DIMENSIONLESS VERSION
//math constants
const pi = 3.141592653589793;
const Msun = 4 * pi**2; //AU
const mau = 6.68458e-12; //meters to AU
const kmkpc = 3.24078e-17;
const kpcau = 2.063 * (10**8); //kpc to au
const aukpc = 4.84814e-9;
const yrsec: real = 60*60*24*365.24; //seconds to years
const secyr: real = 1 / yrsec;
const kpckm = 1.0 / kmkpc;
const G = 6.67430e-11; //m^3 / (kg sec^2)
const G_code = 6.67430e-11 * ((1e-3)**3) * (yrsec)**2 * kmkpc**3 * (1/toCodeTime(1))**2  * toCodeLength(1)**3 * (1/toCodeMass(1)); //G in km^3 / (kg * sec^2)
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
var pot = 6; //set galactic potential to that of imported soliton
var period: real = 0.0;
/* config const po = 27492.260803351877;
config const rc = 0.05359134269304645; */
config const po = (0.5**4) * 27492.260803351877;
config const rc = 2.0*0.05359134269304645;


var Dom = {1..0,1..0,1..0};
//var Dom = {{1..0,1..0,1..0},{1..0,1..0,1..0},{1..0,1..0,1..0}};
var arr : [Dom] real;

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
    acc = force((r0,0.0,0.0),pot,calcpar);
    calcpar[5] = par[1] / r0; //Rhalo becomes dimensionless
    period = 2 * pi * r0 / sqrt(len(acc) * r0); //in seconds; 2pi * r/v

    //period = 2 * pi * r0 / magVel; //in seconds; 2pi * r/v
    calcpar[0] = calcpar[0] * (period ** 2)/ (r0 ** 3); //GM becomes dimensionless
    dt = dt / period; //divide dt by period to make it dimensionless
    writeln("period ", period);
    writeln("dt ", dt);

    }
  else if pot == 0 {
    calcpar[0] = G_code * 1000000000.0 * (1.989e30); //GM in km^3 / sec^2
  }
  else if pot == 6 {
    /* var Dom = {1..0,1..0,1..0};
    var arr : [Dom] real; */
    var ff = openH5File("test/saves/phi_000000.h5");
    readRealArray(ff,'phi',Dom,arr);
    writeln(Dom);

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
  else if pot == 5 || pot == 6 {
    r = len(pos);
    var height: real = 0.0;
    var dr: real = 0.1;
    var mass: real = 0.0;
    var bound: int = (r/dr): int;
    var a: real = 0.0;
    for i in 0..bound {
      const r1 = i * dr;
      height = po * 4 * pi * (r1**2) / (1 + 0.091 * (i * dr/rc)**2)**8; //
      //height = po / (1 + 0.091 * (i * dr/rc)**2)**8; //
      mass = mass + (height * dr);
      //a += -G * height * dr / (i*dr)**2;
    }
    //writeln(mass);
    acc = (- G_code * mass / r**2) * (pos/r);
    //acc = a * (pos/r);
  }
  return acc;
}

//initialize particles in a ring
proc init_ring (pos, vel, AM, SD, PS, calcpar,WritingChannel,AMWritingChannel,PSWritingChannel) {
  //assign counterclockwise from rightmost point on x axis.
  var posAngle, velAngle, magVel: real;
  var acc: 3*real;
  /* acc = force((r0,0.0,0.0),pot,calcpar); */
  //acc = force((1.0,0.0,0.0),pot,calcpar);
  acc = search_force((r0,0.0,0.0));
  //magVel = sqrt(len(acc) * r0) * period / r0;
  magVel = sqrt(len(acc) * r0); //already dimensionless
  //magVel = sqrt(len(acc));

  //var velAngle: real;
  for j in 1..Ne { //Ne=number of desired particles
    posAngle = j * 2.0 * pi / Ne; //cycle thru angles in a circle

    pos[j][0] = r0 * cos(posAngle);
    pos[j][1] = r0 * sin(posAngle);
    pos[j][2] = 0.0;
    /* pos[j][0] = cos(posAngle); //dimensionless
    pos[j][1] = sin(posAngle);
    pos[j][2] = 0.0; */

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

  //update position of jth particle
  pos += dt * vel;
  //writeln(vel);
  //search acceleration
  var a: 3*real;
  //writeln("searching for force");
  a = search_force(pos); //returns dimensionless acc
  //writeln("found force");

  //writeln(a);
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

proc search_force(pos){
  //search acceleration
  var acc: 3*real;
  if pot == 6 {
    var dr: real = 128 / 2.0;
    var x_pos: int = floor(abs(-0.9921875 - pos[0]) * dr): int;
    var y_pos: int = floor(abs(-1.00585938 - pos[1]) * dr): int;
    var z_pos: int = floor(abs(-1.0 - pos[2]) * dr): int;
    //writeln("x pos ",abs(-1.0 - pos[0]) * dr-0.25," y pos ",abs(-1.0 - pos[1]) * dr-0.25," z pos ",abs(-1.0 - pos[2]) * dr-0.25);
    /* writeln("x pos ",x_pos," y pos ",y_pos," z pos ",z_pos); */
    //var a: [0..2] real;
    var a: real;
    a = arr[x_pos,y_pos,z_pos];
    /* writeln("a ",a); */
    //take the derivative of a
    var h: real = 2.0 * (2.0 / 128);

    //derivative blows up at edges, check for edge
    if x_pos > 126 {
      x_pos = 126;
    }
    if y_pos > 126 {
      y_pos = 126;
    }
    if z_pos > 126 {
      z_pos = 126;
    }

    acc(0) = -1*(arr[x_pos + 1, y_pos,z_pos] - arr[x_pos-1,y_pos,z_pos])/ h;
    acc(1) = -1*(arr[x_pos, y_pos + 1,z_pos] - arr[x_pos,y_pos-1,z_pos])/ h;
    acc(2) = -1*(arr[x_pos, y_pos,z_pos + 1] - arr[x_pos,y_pos,z_pos-1])/ h;
    /* writeln("f(x+1) ",arr[x_pos+1, y_pos,z_pos]);
    writeln("f(x-1) ",arr[x_pos-1, y_pos,z_pos]);
    writeln("f(x+1) - f(x-1) ",(arr[x_pos+1, y_pos,z_pos] - arr[x_pos-1, y_pos,z_pos])/(2*0.015625)); */


    //acc = acc * pos /len_pos;

    /* acc(0) = a[0];
    acc(1) = a[1];
    acc(2) = a[2]; */
  }
  else {
    //searches up potential from box and converts to dimensionless units
    var nfwbox = open("nfwbox.csv",iomode.r); // open box
    var boxchannel = nfwbox.reader(); //open reading channel
    var dr: real = 128.0 / toCodeLength(40); //nsteps / range(made dimensionless)

    //writeln(dr);
    var tmp: string;
    //determine column number from z position
    var col: int = floor(abs(toCodeLength(-20.0) - pos[2]) * dr): int;
    var row: real = floor(abs(toCodeLength(-20.0) - pos[0]) * dr);
    //determine row number from x and y positions
    row = 128 * row; //row = 10 * row;
    row += (floor(abs(toCodeLength(-20.0) - pos[1]) * dr));
    var frow: int = row: int; //final row number
    //writeln(frow);
    //writeln(col);
    //writeln(pos);
    for i in 1..frow{
      boxchannel.readln(tmp);
    }
    for i in 0..col{
      boxchannel.read(tmp);
    }
    //writeln(tmp);
    var a: real = tmp: real;
    //acc = -(pos/len(pos)) * a * (period**2)/r0; //acc made dimensionless
    acc = -(pos/len(pos)) * a;
    //writeln(acc);
  }
  return acc;

}
//shifts velocity by a halfstep in direction of sign
proc halfstep(p, ref v, pot, dt, sign, calcpar) {
  var signed_dt: real = sign * dt;
  var a: 3*real;
  a = search_force(p);
  //writeln("dot force and vel ",a[0]*v[0] + a[1]*v[1] + a[2]*v[2]);

  //writeln("signed dt ",signed_dt," calc from halfstep: ",0.5 * signed_dt * a);
  //writeln("a from halfstep ",a);
  v += 0.5 * signed_dt * a;
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
//update positions at every timestep based on potential box

proc toSeconds (l) { //not in use
  var a: real = ((8 * pi)/(3* H0**2 *omegaM0))**0.5;
  return a / l**2;
}

proc toMeters (l) { //not in use
  var a: real = (8 * pi * (hbar**2))/(3 * (maxion**2) * (H0**2) * omegaM0);
  return (a**0.25) / l;
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

/*

proc getPot(pos) {
  convert pos to code units
  get potential value from position in code units
  take gradient in code units to convert to acceleration
  convert acceleration in code units to regular units
}
*/



proc main(){

  //create array of tuples to hold position and velocity and velocity at each timestep
  var pos: [1..Ne] 3*real;
  var vel: [1..Ne] 3*real;
  var AM: [1..Ne] real; //stores angular momentum of each particle, updated at each timestep
  var SD: [0..N] real; //stores std dev at each timestep
  var PS: [1..nbins] real; //whole array gets updated at each timestep
  //hardcode galactic potential parameters
  var calcpar: [0..5] real;

  var chfile = open("mathematica/bp19.csv",iomode.cw); //create test.csv and open
  var SDfile = open("mathematica/SD19.csv",iomode.cw); //create test.csv and open
  var AMfile = open("mathematica/AM19.csv",iomode.cw); //create test.csv and open
  var PSfile = open("mathematica/PS19.csv",iomode.cw); //create test.csv and open
  var PSWritingChannel = PSfile.writer(); //open writing channel to test.csv
  var SDWritingChannel = SDfile.writer(); //open writing channel to test.csv
  var WritingChannel = chfile.writer(); //open writing channel to test.csv
  var AMWritingChannel = AMfile.writer(); //open writing channel to test.csv

  load_pot(calcpar,pot);
  //set initial position and velocity of every particle in ring
  init_ring(pos, vel, AM, SD, PS, calcpar, WritingChannel,AMWritingChannel,PSWritingChannel);
  //dt = dt / period;
  writeln("initialized ring");
  writeln("pos 1 ",pos[1]);
  writeln("vel 1", vel[1]);
  //
  /* writeln("calling search force ",search_force((0.5,0.5,0.5))); */
  //writeln("calling search force ",search_force((0.4,0.4,0.4)));

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
