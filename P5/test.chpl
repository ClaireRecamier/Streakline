use Random;
use IO;
var rfile = open("ran3.csv",iomode.cw); //create test.csv and open
var channel = rfile.writer();
var randStream = new RandomStream(real);
var randStreamSeeded = new RandomStream(real, 0);
var randsFromStream: [1..1000] real;
randStream.fillRandom(randsFromStream);
for i in 1..1000 {
  randsFromStream[i] = 2.0 * (randsFromStream[i] - 0.5);
  channel.write(randsFromStream[i],"\n");
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
