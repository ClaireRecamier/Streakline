use Random;
use IO;
var rfile = open("nran2.csv",iomode.cw); //create test.csv and open
var channel = rfile.writer();
var randStream = new RandomStream(real);
var randStreamSeeded = new RandomStream(real, 0);
var randsFromStream: [1..2] real;

var r1, r2: 3*real;
nran(r1,r2);

for i in 1..5000 {
  nran(r1,r2);
  channel.write(r1(0),"\n",r1(1),"\n",r1(2),"\n");
  channel.write(r2(0),"\n",r2(1),"\n",r2(2),"\n");

}

proc nran (ref r1, ref r2) {

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


/*
for i in 1..10000 {
  if randsFromStream[i] < 0 then writeln(randsFromStream[i]);
}

randStreamSeeded.fillRandom(randsFromStream);
for i in 1..10000 {
  if randsFromStream[i] < 0 then writeln(randsFromStream[i]);
}
*/
