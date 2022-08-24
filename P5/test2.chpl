use IO;
var cfile2 = open("test7.csv",iomode.r); //open position offsets
var myReadingChannelp = cfile2.reader(); //open reading channel to test.csv
var tmp: string;
myReadingChannelp.read(tmp);
writeln(tmp);
myReadingChannelp.read(tmp);
writeln(tmp);
