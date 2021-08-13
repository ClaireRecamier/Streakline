use IO;
use Math;
use Random;

const pi = 3.141592653589793;
const Msun = 4 * pi**2; //AU
const mau = 6.68458e-12; //meters to AU
const kpcau = 2.063 * (10**8); //kpc to au
const aukpc = 4.84814e-9;
const secyr: real = 60*60*24*365.24; //seconds to years
config const integrator = 0;   //set Integrator to LF
config const N = 6000;//6000 set number of timesteps to equal total of 6 Gyr
config const dt = 1000000.0; //set timestep to Myr
config const mcli: real = 20000 * Msun; //initial mass of cluster
config const mclf: real = 20000 * Msun; //final mass of cluster
config const M = 1; //particles are released every Mth timestep
config const Rcl = 20 * 0.001 * kpcau; //radius of plummer core
var pot = 3; //set galactic potential to that of NFW
var Ne = ceil(N/M) : int; //number of particles released
var k: int = 0; //record how many particles have been released
var dm = (mcli - mclf)/N; //amount of mass released per timesteps
var mcl: real = mcli; //current mass of cluster
var chfile = open("test24.csv",iomode.cw); //create test.csv and open
var myWritingChannel = chfile.writer(); //open writing channel to test.csv
var offset: [0..1] real = [0.5,0.5];
var calcpar: [0..5] real;
var randStream = new RandomStream(real);

proc main () {
  //create array of tuples to hold position and velocity and velocity at each timestep
  var pos: [0..N] 3*real;
  var vel: [0..N] 3*real;

  var pos_lead: [1..Ne] 3*real;
  var pos_trail: [1..Ne] 3*real;
  var vel_lead: [1..Ne] 3*real;
  var vel_trail: [1..Ne] 3*real;

  //hardcode galactic potential parameters
  var par: [0..5] real = [417.0 * (10**3) * mau * secyr, 36.54 * kpcau, 90.0 * pi / 180, 1.0, 1.0 ,0.94];
  if pot == 3 { //if using triaxial NFW potential
    //assuming par = [V, rhalo, phi, q_1, q_2, q_z]
    //calcpar = [GM, c1, c2, c3, c4, rhalo]
    var cosphi: real = cos(par[2]);
    var sinphi: real = sin(par[2]);
    calcpar[0] = par[0]*par[0]*par[1]; //GM
    calcpar[1] = (cosphi**2)/(par[3]*par[3]) + (sinphi**2)/(par[4]*par[4]);
    calcpar[2] = (cosphi**2)/(par[4]*par[4]) + (sinphi**2)/(par[3]*par[3]);
    calcpar[3] = 2*sinphi*cosphi*(1/(par[3]**2) - 1/(par[4]**2));
    calcpar[4] = 1/(par[5]*par[5]);
    calcpar[5] = par[1];
    //writeln("initializing param ",calcpar);
  }


  //hardcode initial position and velocity
  pos[0]=(50.0*kpcau,0,0);
  vel[0]=(0,0.5* sqrt(calcpar[0]/len(pos[0])),0); //(if pot 3) set velocity equal to half of centripetal velocity
  //vel[0]=(0,0.5* sqrt(1000000000*Msun/len(pos[0])),0); // (if pot 0)set velocity equal to half of centripetal velocity

  myWritingChannel.write("x cluster,y cluster,x cluster vel,y cluster vel,x lead trail,y lead tail,x vel lead tail, y vel lead tail,x trail tail,y trail tail,x vel trail tail, y vel trail tail\n");

  //writeln("initial energy ", energy(pos,vel,0));
  back_orbit(pos, vel, pot, integrator, N, dt, pos_lead, pos_trail, vel_lead, vel_trail);
  pos[0] = pos[N-1];
  vel[0] = vel[N-1];
  fwd_orbit(pos, vel, pot, integrator, N, dt, pos_lead, pos_trail, vel_lead, vel_trail);
  //writeln("final energy ", energy(pos,vel,N));
  myWritingChannel.close();
  chfile.close();
}

//orbit procedure: advances cluster in position and velocity using integrator of choice by N timesteps
proc fwd_orbit (pos, vel, pot, integrator, N, dt, pos_lead, pos_trail, vel_lead, vel_trail) {
    if integrator == 0  { //if leapfrog
      //move velocity forward half a timestep
      halfstep(pos[0],vel[0],pot,dt,1.0);
      //myWritingChannel.write(aukpc * pos[0][0],",",aukpc * pos[0][1],",",aukpc * vel[0][0],",",aukpc * vel[0][1],"\n");
      //writeln("pos cluster after first halfstep: ",pos[0]);
      //writeln("vel cluster after first halfstep: ",vel[0]);
      for i in 1..N-1 {//make N full steps in pos and vel forwards
        //mcl -= dm;//decrease mass

        leapfrog(pos,vel,i,dt,1.0);
        //myWritingChannel.write(aukpc * pos[i][0],",",aukpc * pos[i][1],",",aukpc * vel[i][0],",",aukpc * vel[i][1],",");

        for j in 1..k {
          stream_step(pos_lead[j], vel_lead[j], pos[i], dt);
          stream_step(pos_trail[j], vel_trail[j], pos[i], dt);
          //myWritingChannel.write(",",pos_lead[j][0],",",pos_lead[j][1],",",vel_lead[j][0],",",vel_lead[j][1],",",pos_trail[j][0],",",pos_trail[j][1],",",vel_trail[j][0],",",vel_trail[j][1]);
          //writeln("stars ",j," at ",i," step ",aukpc * pos_lead[j][0],",",aukpc * pos_lead[j][1],",",aukpc * vel_lead[j][0],",",aukpc * vel_lead[j][1],",",aukpc * pos_trail[j][0],",",aukpc * pos_trail[j][1],",",aukpc * vel_trail[j][0],",",aukpc * vel_trail[j][1]);
        }

        if i % M == 0 {
          k+=1;
          eject(pos[i],vel[i],pos_lead[k], vel_lead[k], pos_trail[k], vel_trail[k]);
          //myWritingChannel.write(aukpc * pos_lead[k][0],",",aukpc * pos_lead[k][1],",",aukpc * vel_lead[k][0],",",aukpc * vel_lead[k][1],",",aukpc * pos_trail[k][0],",",aukpc * pos_trail[k][1],",",aukpc * vel_trail[k][0],",",aukpc * vel_trail[k][1],"\n");
          //writeln("after ejecting ", vel_lead[k]);
          //writeln(pos_trail[k]);
        }
        /* //print snapshot at 5000th timestep
        if i == N-1000 {
          for j in 1..k {
            myWritingChannel.write(pos_lead[j][0],",",pos_lead[j][1],",",vel_lead[j][0],",",vel_lead[j][1],",",pos_trail[j][0],",",pos_trail[j][1],",",vel_trail[j][0],",",vel_trail[j][1],"\n");
          }
        }
        */
        //myWritingChannel.write("\n");

      }
      //move velocity backward half a timestep
      halfstep(pos[N-1],vel[N-1],pot,dt,-1.0);
      //position of cluster at last timestep
      //myWritingChannel.write(aukpc * pos[N-1][0],",",aukpc * pos[N-1][1],",",aukpc * vel[N-1][0],",",aukpc * vel[N-1][1],"\n");

      //positions of streams at last timestep
      for j in 1..k {
        myWritingChannel.write(aukpc * pos[j][0],",",aukpc * pos[j][1],",",aukpc * vel[j][0],",",aukpc * vel[j][1],",");
        myWritingChannel.write(aukpc * pos_lead[j][0],",",aukpc * pos_lead[j][1],",",aukpc * vel_lead[j][0],",",aukpc * vel_lead[j][1],",",aukpc * pos_trail[j][0],",",aukpc * pos_trail[j][1],",",aukpc * vel_trail[j][0],",",aukpc * vel_trail[j][1],"\n");
      }
  }

}

proc back_orbit (pos, vel, pot, integrator, N, dt, pos_lead, pos_trail, vel_lead, vel_trail)
{
  if integrator == 0  { //if leapfrog
    //move velocity forward half a timestep
    //writeln("pos cluster before first halfstep: ",aukpc * pos[0]);
    //writeln("vel cluster before first halfstep: ",aukpc * vel[0]);
    halfstep(pos[0],vel[0],pot,dt,-1.0);
    //myWritingChannel.write(aukpc * pos[0][0],",",aukpc * pos[0][1],",",aukpc * vel[0][0],",",aukpc * vel[0][1],"\n");
    //writeln("pos cluster after first halfstep: ",aukpc * pos[0]);
    //writeln("vel cluster after first halfstep: ",aukpc * vel[0]);
  for i in 1..N-1 {//make N-1 full steps in pos and vel forwards
    //mcl -= dm;//decrease mass

    leapfrog(pos,vel,i,dt,-1.0);
    //myWritingChannel.write(aukpc * pos[i][0],",",aukpc * pos[i][1],",",aukpc * vel[i][0],",",aukpc * vel[i][1],"\n");
    //writeln(aukpc * pos[i][0],",",aukpc * pos[i][1],",",aukpc * vel[i][0],",",aukpc * vel[i][1]);
  }
  //move velocity backward half a timestep
  halfstep(pos[N-1],vel[N-1],pot,dt,1.0);
  //myWritingChannel.write(aukpc * pos[N-1][0],",",aukpc * pos[N-1][1],",",aukpc * vel[N-1][0],",",aukpc * vel[N-1][1],"\n");
  //writeln("ending pos ",pos[N-1]," ending vel ",vel[N-1]);
  }
}

//procedure to advance ejected particles by a timestep
proc stream_step(ref pos, ref vel, pos_cl, dt) {
  //update position of jth particle
  pos += dt * vel;
  //calculate acceleration from Msun
  var a: 3*real;
  a = force(pos,pot);
  //writeln("pos ",pos);
  //writeln("a ",a);
  //writeln("force plummer ",force_plummer(pos_cl - pos, Rcl));
  //calculate plummer acceleration
  a += force_plummer(pos_cl - pos, Rcl);
  //update velocity of jth particle
  vel += dt * a;
}

//procedure to eject particles
proc eject(pos_cl, vel_cl, ref pos_lead, ref vel_lead, ref pos_trail, ref vel_trail)
{

  //initialize velocity offsets
  var dvl, dvt: real;

  var ran_vec: 3*real = norm_rand(); //get new tuple of normalized random numbers
  dvl = len(ran_vec)*offset[1]/3; //convert to maxwell distribution

  ran_vec = norm_rand(); //get new tuple of normalized randomn\ numbers
  dvt = len(ran_vec)*offset[1]/3; //convert to maxwell distribution
  //myWritingChannel.write(dvl,",",dvt,",");

  //calculate angular velocity of Cluster
  //writeln("pos cl at ejection ",pos_cl, " vel cl at ejection", vel_cl);
  var omega: 3*real;
  omega = cross(pos_cl, vel_cl); //r x v
  //writeln("omega of cluster: ",omega);
  var r = len(pos_cl);
  omega = omega / (r**2); //(r x v)/r^2
  var om: real = len(omega);

  var Rj: 3*real = tidal_radius(pos_cl, vel_cl, om);//calculate tidal radius
  var dRRj: real = len(Rj) * offset[0];
  //initial position of particle is position of cluster plus or minus tidal radius
  //writeln("Rj: ",Rj);
  //pos_lead = pos_cl - (Rj * (pos_cl/r)); //multiplied by unit vector
  //pos_trail = pos_cl + (Rj * (pos_cl/r));


  ran_vec = norm_rand();//get new tuple of normalized random numbers
  pos_lead = pos_cl - Rj + (dRRj * ran_vec);

  ran_vec = norm_rand(); //get new tuple of normalized random numbers
  pos_trail = pos_cl + Rj + (dRRj * ran_vec);

  //pos_lead = pos_cl - Rj;
  //pos_trail = pos_cl + Rj;

  var mag_vcl: real = len(vel_cl);
  var uv_vel: 3*real = vel_cl / mag_vcl; // unit vector of vel vector
  var uv_pos: 3*real = pos_cl / r; //unit vector of pos vector

  vel_lead = uv_vel * (mag_vcl - len(Rj) * om) - (dvl * uv_pos);
  vel_trail = uv_vel * (mag_vcl + len(Rj) * om) + (dvt * uv_pos);
  /*
  var mag_vcl: real = len(vel_cl);
  var uv_vel: 3*real = vel_cl / mag_vcl; // unit vector of vel vector
  vel_lead = uv_vel * (mag_vcl - len(Rj) * om);
  vel_trail = uv_vel * (mag_vcl + len(Rj) * om);
  */

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
  //writeln("signed dt ",signed_dt," calc from halfstep: ",0.5 * signed_dt * a);
  //writeln("a from halfstep ",a);
  v += 0.5 * signed_dt * a;
}

//shifts velocity and position by a full step forward
proc leapfrog(pos, vel, i, dt, sign) {
  var a: 3*real;
  pos[i] = pos[i-1] + (vel[i-1] * dt * sign);
  a = force(pos[i],pot);
  vel[i] = vel[i-1] + (dt * sign * a);
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
  var r: real;
  if pot == 0 { //if using point mass potential
    r = len(pos);
    //acc = (-1,-1,-1)*(Msun * pos)/dist**3; //assumes the sun stays at origin
    acc = (-1,-1,-1)*(100000000000 * Msun * pos)/r**3;
    //writeln("acc ",acc);
  }
  else if pot == 3 { //NFW triaxial potential, C version
    //calcpar = [GM, c1, c2, c3, c4, rhalo]
    r = sqrt(calcpar[1]*pos[0]*pos[0] + calcpar[2]*pos[1]*pos[1] + calcpar[3]*pos[0]*pos[1] + calcpar[4]*pos[2]*pos[2]);
    var aux: real = 0.5 * calcpar[0] / (r**3) * (1.0/(1.0 + calcpar[5]/r)-log(1.0+r/calcpar[5]));

    acc[0]=aux*(2*calcpar[1]*pos[0] + calcpar[3]*pos[1]);
    acc[1]=aux*(2*calcpar[2]*pos[1] + calcpar[3]*pos[0]);
    acc[2]=aux*(2*calcpar[4]*pos[2]);
    //writeln("acc ",acc);
  }
  return acc;
}

proc force_plummer(r, Rcl) {
  // dist from cluster to particle, cluster radius
  var dist = len(r);
  var raux: real = sqrt(dist**2 + Rcl**2);
  var acc: 3*real;
  acc = (mcl / raux**3 ) * r;
  //write("{",acc[0],",",acc[1],"},");
  return acc;
}

proc tidal_radius(pos_cl,vel_cl,omega) {
  //writeln("pos cl ",pos_cl," vel cl ",vel_cl," omega ",omega);
  var un_vec: 3*real = pos_cl/len(pos_cl);
  var delta = 0.02 * kpcau * un_vec;
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
  //writeln("mcl ",mcl);
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

proc norm_rand () {
  var arr: [0..2] real;
  var rands: 3*real;
  randStream.fillRandom(arr);
  for i in 0..2 {
    rands[i] = 3.5 * (arr[i] - 0.5);
  }
  return rands;
}
