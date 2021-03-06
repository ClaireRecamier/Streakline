use IO;
use Math;
use Random;
//DIMENSIONFULL version
//math constants
const pi = 3.141592653589793;
const Msun = 4 * pi**2; //AU
const mau = 6.68458e-12; //meters to AU
const kpcau = 2.063 * (10**8); //kpc to au
const aukpc = 4.84814e-9;
const secyr: real = 60*60*24*365.24; //seconds to years
//setup constants
config const integrator = 0;   //set Integrator to LF
config const N = 10;//6000 set number of timesteps to equal total of 6 Gyr
config const dt = 1000000.0; //set timestep to Myr
config const offsetOption = 0; //0 = no radial offsets, 1 = using C generated random numbers, 2 = generating own random numbers
config const r0 = 15 * kpcau; // radius of ring, in kpc
config const Ne = 1; //number of particles in ring
config const nbins = 4;
var pot = 3; //set galactic potential to that of NFW

var k: int = 0; //record how many particles have been released
var chfile = open("test8.csv",iomode.cw); //create test.csv and open
var SDfile = open("SD8.csv",iomode.cw); //create test.csv and open
var AMfile = open("AM8.csv",iomode.cw); //create test.csv and open
var PSfile = open("PS8.csv",iomode.cw); //create test.csv and open
var PSWritingChannel = PSfile.writer(); //open writing channel to test.csv
var SDWritingChannel = SDfile.writer(); //open writing channel to test.csv
var WritingChannel = chfile.writer(); //open writing channel to test.csv
var AMWritingChannel = AMfile.writer(); //open writing channel to test.csv
var offset: [0..1] real = [0.2,0.2];
var randStream = new RandomStream(real);

proc main () {
  //create array of tuples to hold position and velocity and velocity at each timestep
  var pos: [1..Ne] 3*real;
  var vel: [1..Ne] 3*real;
  var AM: [1..Ne] real;
  var SD: [0..N] real;
  var PS: [1..nbins] real;

  //hardcode galactic potential parameters
  var calcpar: [0..5] real;
  var par: [0..5] real = [417.0 * (10**3) * mau * secyr, 36.54 * kpcau, 90.0 * pi / 180, 1.0, 1.0 ,1.0];
  if pot == 3 { //if using triaxial NFW potential
    //assuming par = [V, rhalo, phi, q_1, q_2, q_z]
    //calcpar = [GM, c1, c2, c3, c4, rhalo]
    var cosphi: real = cos(par[2]);
    var sinphi: real = sin(par[2]);
    calcpar[0] = par[0]*par[0]*par[1]; //GM
    //writeln(calcpar[0]);
    calcpar[1] = (cosphi**2)/(par[3]*par[3]) + (sinphi**2)/(par[4]*par[4]);
    calcpar[2] = (cosphi**2)/(par[4]*par[4]) + (sinphi**2)/(par[3]*par[3]);
    calcpar[3] = 2*sinphi*cosphi*(1/(par[3]**2) - 1/(par[4]**2));
    calcpar[4] = 1/(par[5]*par[5]);
    calcpar[5] = par[1];
    //writeln("initializing param ",calcpar);
    //writeln("c1 " , calcpar[1], " c2 ",calcpar[2]," c3 ",calcpar[3], " c4 ",calcpar[4]);
  }

  //initialize initial position and velocity of every particle in ring
  init_ring(pos, vel, AM, SD, PS, calcpar);

  fwd_orbit(pos, vel, AM, SD, PS, pot, integrator, N, dt, calcpar);

  WritingChannel.close();
  SDWritingChannel.close();
  AMWritingChannel.close();
  PSWritingChannel.close();
  PSfile.close();
  AMfile.close();
  SDfile.close();
  chfile.close();
}

proc init_ring (pos, vel, AM, SD, PS, calcpar) {
  //assign counterclockwise from rightmost point on x axis.
  //3000 particles in bottom and 30000 on top
  var posAngle, velAngle, magVel: real;
  //var velAngle: real;
  for j in 1..Ne {
    posAngle = j * 2.0 * pi / Ne;
    pos[j][0] = r0 * cos(posAngle);
    pos[j][1] = r0 * sin(posAngle);
    pos[j][2] = 0.0;
    velAngle = (pi / 2.0) - posAngle;
    magVel =  sqrt(calcpar[0] / r0);
    //magVel = sqrt(100000000000*Msun/r0);
    vel[j][0] = magVel * cos(velAngle);
    vel[j][1] = -1.0 * magVel * sin(velAngle);
    vel[j][2] = 0.0;

    AM[j] = angMom(pos[j],vel[j]);
    writeln("initial pos ", pos[j] * aukpc);
    writeln("initial vel ", vel[j] * aukpc / secyr);
    WritingChannel.write((pos[j][0]*vel[j][0] + pos[j][1]*vel[j][1] + pos[j][2]*vel[j][2]) * aukpc**2,",0.0,0.0,",len(vel[j])*aukpc,",",len(pos[j])*aukpc,",");
    WritingChannel.write(pos[j][0] * aukpc,",",pos[j][1] * aukpc,",",pos[j][2] * aukpc,",");
    WritingChannel.write(vel[j][0] * aukpc ,",",vel[j][1] * aukpc,",",vel[j][2] * aukpc,"\n");
    AMWritingChannel.write(AM[j],",");

  }
  pow_spec(pos,vel,PS);
  AMWritingChannel.write("\n");

  //iterate through 2pi radians divided by 6000 particles.
  //for each, multiply r0 by sin to get y, r0 by cosine to get x
  //calculate magnitude of centripetal velocity vector for each particles
  //angle velocity vector makes with axis is 90 - theta.
  //multiple magnitude of velocity vector with sin and cosine of angle.
}

proc pow_spec (pos,vel, PS) {//gets computed at each timstep
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

//orbit procedure: advances cluster in position and velocity using integrator of choice by N timesteps
proc fwd_orbit (pos, vel, AM, SD, PS, pot, integrator, N, dt, calcpar) {
    if integrator == 0  { //if leapfrog
      //move velocity of all particles a halfstep forward
      /*
      for j in 1..Ne {
        halfstep(pos[j],vel[j],pot,dt,1.0, calcpar);
        AM[j] = angMom(pos[j],vel[j]);
        //writeln("vel after halfstep ",vel[j]);
        //writeln("pos after halfstep ",pos[j]);


        //WritingChannel.write(pos[j][0] * aukpc,",",pos[j][1] * aukpc,",",pos[j][2] * aukpc,",");
        //WritingChannel.write(vel[j][0] * aukpc,",",vel[j][1] * aukpc,",",vel[j][2] * aukpc,"\n");
        //AMWritingChannel.write(AM[j],",");
      }
      */


      //AMWritingChannel.write("\n");
      for i in 0..N { //for each subsequent timestep
        //move velocity of all particles a fullstep forwards
        for j in 1..Ne {
          stream_step(pos[j], vel[j], dt, calcpar);
          AM[j] = angMom(pos[j],vel[j]);
          //writing position and velocity of to file (one particle per line)
          WritingChannel.write(pos[j][0] * aukpc,",",pos[j][1] * aukpc,",",pos[j][2] * aukpc,",");
          WritingChannel.write(vel[j][0] * aukpc,",",vel[j][1] * aukpc,",",vel[j][2] * aukpc,"\n");
          //writing angular momentum to a diff file
          AMWritingChannel.write(AM[j],",");
          }
        AMWritingChannel.write("\n");
        SD[i] = stdev(AM);
        pow_spec(pos,vel,PS);
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
        SDWritingChannel.write(i,",",SD[i],"\n");
      }
  }
}

proc angMom(pos,vel){
  var omega: 3*real;
  omega = cross(pos, vel); //r x v
  omega = omega / (len(pos)**2); //(r x v)/r^2
  return len(omega);
}

//procedure to advance ejected particles by a timestep. only subject to force from central pot
proc stream_step(ref pos, ref vel, dt, calcpar) {
  //update position of jth particle
  //writeln("dot pos and vel before ",pos[0]*vel[0] + pos[1]*vel[1] + pos[2]*vel[2]);
  pos += dt * vel;
  //calculate acceleration from Msun
  var a: 3*real;
  a = force(pos,pot, calcpar);


  //writeln("force and pos parallel? ",cross(pos,a));
  //update velocity of jth particle
  vel += dt * a;
  //writeln("position ",pos);
  //writeln("velocity ", vel);
  WritingChannel.write(pos[0]*vel[0] + pos[1]*vel[1] + pos[2]*vel[2],",");
  WritingChannel.write(a[0]*vel[0] + a[1]*vel[1] + a[2]*vel[2],",");
  WritingChannel.write(len(a),",");
  WritingChannel.write(len(vel),",");
  WritingChannel.write(len(pos),",");

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
  return SD/r0;
}
//shifts velocity by a halfstep in direction of sign
proc halfstep(p, ref v, pot, dt, sign, calcpar) {
  var signed_dt: real = sign * dt;
  var a: 3*real;
  a = force(p,pot,calcpar);
  writeln("dot force and vel ",a[0]*v[0] + a[1]*v[1] + a[2]*v[2]);

  //writeln("signed dt ",signed_dt," calc from halfstep: ",0.5 * signed_dt * a);
  //writeln("a from halfstep ",a);
  v += 0.5 * signed_dt * a;
}


proc force(pos,pot,calcpar){
  var acc: 3*real;
  var r: real;
  if pot == 0 { //if using point mass potential
    //r = len(pos);
    //acc = (-1,-1,-1)*(Msun * pos)/dist**3; //assumes the sun stays at origin
    //acc = (-1,-1,-1)*(100000000000 * Msun * pos)/r**3;
    //writeln("acc ",acc);
    acc = (-1,-1,-1)*(100000000000 * Msun * pos)/r0**3;

  }
  else if pot == 3 { //NFW triaxial potential, C version
    //calcpar = [GM, c1, c2, c3, c4, rhalo]
    r = sqrt(calcpar[1]*pos[0]*pos[0] + calcpar[2]*pos[1]*pos[1] + calcpar[3]*pos[0]*pos[1] + calcpar[4]*pos[2]*pos[2]);
    var aux: real = 0.5 * calcpar[0] / (r**3) * (1.0/(1.0 + calcpar[5]/r)-log(1.0+r/calcpar[5]));

    acc[0]=aux*(2*calcpar[1]*pos[0] + calcpar[3]*pos[1]);
    acc[1]=aux*(2*calcpar[2]*pos[1] + calcpar[3]*pos[0]);
    acc[2]=aux*(2*calcpar[4]*pos[2]);
    writeln("acc magnitude ", len(acc));
  }
  return acc;
}


proc len((x,y,z)) {
  return sqrt(x**2 + y**2 + z**2);
}

proc toreal (str) {
  //writeln(str);
  var numb: [0..1] int;
  var n = str.split(".");
  numb[0] = abs(n[0]:int);
  numb[1] = n[1]:int;
  var n1: real = 1.0 * numb[1] / 10**n[1].size + numb[0];
  if n[0][0] == '-' then n1 = n1 * -1.0;
  //writeln(n1);
  return n1;
}

proc norm_rand (ref r1, ref r2, randStream) {
  //fills two tuples with normally distributed random numbers
  var arr1: [0..2] real;
  var arr2: [0..2] real;
  randStream.fillRandom(arr1);
  randStream.fillRandom(arr2);

  r1[0] = sqrt(-2 * log(arr1[0])) * cos(2*pi*arr1[1]);
  r1[1] = sqrt(-2 * log(arr1[0])) * sin(2*pi*arr1[1]);

  r1[2] = sqrt(-2 * log(arr1[2])) * cos(2*pi*arr2[0]);
  r2[0] = sqrt(-2 * log(arr1[2])) * sin(2*pi*arr2[0]);

  r2[1] = sqrt(-2 * log(arr2[1])) * cos(2*pi*arr2[2]);
  r2[2] = sqrt(-2 * log(arr2[1])) * sin(2*pi*arr2[2]);
}

proc readTuple(ReadingChannel, ref r){
  //fills tuple r with numbers from reading channel
  var tmp: string;
  for i in 0..2 {
    ReadingChannel.read(tmp);
    r(i) = toreal(tmp);
  }
}

proc loadOffsets(ref dvl, ref dvt, ref r1, ref r2, ReadingChannelv, ReadingChannelp, randStream) {
  if offsetOption == 0 { //if adding no random radial offset
    dvl = 0;
    dvt = 0;
    r1 = (0.0,0.0,0.0);
    r2 = (0.0,0.0,0.0);
    offset = [0.0,0.0];
  }
  else if offsetOption == 1 { //if using C's random radial offsets
    readTuple(ReadingChannelv, r1);
    dvl = len(r1) * offset[1]/3;
    readTuple(ReadingChannelv,r2);
    dvt = len(r2) * offset[1]/3;
    readTuple(ReadingChannelp,r1);
    readTuple(ReadingChannelp,r2);
  }
  else if offsetOption == 2 { //if adding chapel generated random radial offsets
    norm_rand(r1,r2,randStream); //get new tuple of normalized random numbers
    dvl = len(r1)*offset[1]/3; //convert to maxwell distribution
    dvt = len(r2)*offset[1]/3; //convert to maxwell distribution
    norm_rand(r1, r2, randStream);//get new tuple of normalized random numbers
  }
  WritingChannel.write(dvl,",",dvt,",");
}

/*
compute power spectra
keep track of number of particles at all timestep
average is the number of particles divided by number of bins

if you do a noncircular orbit but spherical potential and without radial OFFSETS
 do you still get a cold thin stream?

tRy example woithout a globular cluster, just stream particles in a ring.
compare with pointmass potential and NFW triaxial. how does adding in triaxiality
change how they move? */
/*
how to compute velocity dispersions?
Velocity dispersion is standard deviation of velocities. This is
calculated by dividing standard deviation of angular momentum
For a given point in time:
calculate angular momentum of each particle in stream
find standard deviation of angular momentum
divide standard deviation of angular momentum by r0
*/
