use Random;
use IO;
var rfile = open("RadialOffsets/nran3.csv",iomode.cw); //create test.csv and open
var channel = rfile.writer();
var randStream = new RandomStream(real);
var randStreamSeeded = new RandomStream(real, 0);

var randFromStream: [1..10000] real;
randStreamSeeded.fillRandom(randFromStream);

for i in 1..10000 do channel.write(randFromStream[i],"\n");
