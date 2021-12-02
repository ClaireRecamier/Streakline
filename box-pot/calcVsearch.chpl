use IO;
//UNITS
const pi = 3.141592653589793;
const Msun = 4 * pi**2; //AU
const mau = 6.68458e-12; //meters to AU
const kmkpc = 3.24078e-17;
const kpcau = 2.063 * (10**8); //kpc to au
const aukpc = 4.84814e-9;
const secyr: real = 60*60*24*365.24; //seconds to years
const kpckm = 1.0 / kmkpc;
const G = 6.67430e-11 * ((1e-3)**3); //G in km^3 / (kg * sec^2)

proc pot_box(calcpar,pot) { //creates acceleration box in km/seconds squared
  var nfwbox = open("nfwbox.csv",iomode.cw); //create test.csv and open
  var boxchannel = nfwbox.writer(); //open writing channel to test.csv

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
    }
  else if pot == 0 {
    calcpar[0] = G * 1000000000.0 * (1.989e30); //GM in km^3 / sec^2
  }

  //create box
  //var dr: real = (40.0 * kpckm)/128.0; //in kpc : range / nsteps
  var dr: real = (40.0 * kpckm)/128.0; //in km : range / nsteps
  var x, y, z: real = 0.0;
  var acc: real;
  var row: int = 0;
  var col: int = 0;
  for i in -64..63 {//x
    x = i * dr; //position in km
    for j in -64..63 {//y
      //y = i * dr; //in km
      y = j * dr;
      for k in -64..63 {//z
        z = k * dr; //in km
        //writeln(x,",",y,",",z);
        acc = force((x,y,z),pot,calcpar);
        boxchannel.write(acc," ");
        if i == -39 && j == 0 && k == 0 {
          //writeln(x / (kpckm * 15.0));
          //writeln(y / (kpckm * 15.0));
          //writeln(z / (kpckm * 15.0));
          //writeln(j);
          //writeln("when printing box: ",acc);
          //writeln("row ",row);
          //writeln("col" ,col);
        }
        col += 1;
      }
      col = 0;
      row += 1;
      boxchannel.write("\n");
    }
  }

  boxchannel.close();
  nfwbox.close();
}


proc search_force(pos,calcpar) {
  var nfwbox = open("nfwbox.csv",iomode.r); //create test.csv and open
  var boxchannel = nfwbox.reader(); //open reading channel
  var dr: real = 128.0 / (40.0 / 15.0); //nsteps / range(made dimensionless)
  //writeln(dr);
  var tmp: string;

  writeln("section of z pos ",abs((-20.0/15.0) - pos[2]) * dr);
  writeln("section of x pos ",abs((-20.0/15.0) - pos[0]) * dr);
  writeln("section of y pos ",abs((-20.0/15.0) - pos[1]) * dr);
  //search acceleration
  var acc: 3*real;
  //determine column number from z position
  var col: int = floor(abs((-20.0/15.0) - pos[2]) * dr): int;
  //determine row number from x position
  var row: real = floor(abs((-20.0/15.0) - pos[0]) * dr);

  //find and print input position to target column and corresponding force
  var position: 3*real = (floor(abs((-20.0/15.0) - pos[0]) * dr)*(1.0/dr),floor(abs((-20.0/15.0) - pos[1]) * dr)*(1.0/dr),floor(abs((-20.0/15.0) - pos[2]) * dr)*(1.0/dr));
  writeln((position) - (20.0 /15.0)); //print position in dimensionless units
  writeln(force((position * 15.0 * kpckm) - (20.0 * kpckm),3,calcpar));

  //writeln(abs((-20.0/15.0) - pos[0]) * dr);
  //determine row number from x and y positions
  //writeln(row);
  row = 128 * (row);
  //writeln(row);
  //row += (floor(abs((-20.0/15.0) - pos[1]) * dr));
  row += (floor(abs((-20.0/15.0) - pos[1]) * dr));
  //writeln(row);

  var frow: int = row: int;
  //writeln(frow);
  //writeln(col);
  //writeln(pos);
  //for i in 0..(frow-1){
  for i in 1..frow {
    boxchannel.readln(tmp);
  }
  for i in 0..col {
    boxchannel.read(tmp);
  }
  //writeln(tmp);
  var a: real = tmp: real;
  //acc = (pos/len(pos)) * a * (period**2)/r0; //acc made dimensionless
  //writeln(acc);
  return a;

}

proc len((x,y,z)) {
  return sqrt(x**2 + y**2 + z**2);
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
  return len(acc);
}
proc main(){
  var pot: int = 3;
  var par: [0..5] real = [417.0, 36.54 * kpckm, 90.0 * pi / 180, 1.0, 1.0 ,1.0];
  var calcpar: [0..5] real;
  pot_box(calcpar,pot);
  var pos: 3*real;
  pos = (-0.8,0.0,0.0);
  writeln("force from searching up position in box: ",search_force(pos,calcpar)); //in km/s
  pos = (-0.8 * 15.0 * kpckm, 0.0 * 15.0 * kpckm , 0.0 * 15.0 * kpckm);
  writeln("force from force function: ",force(pos,pot,calcpar));
}
