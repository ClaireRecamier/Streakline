use IO;

var cfile = open("ctest6.csv",iomode.r); //create test.csv and open
var myReadingChannel = cfile.reader(); //open reading channel to test.csv
var r1,r2,r3: string;
myReadingChannel.read(r1);
var arr: string = "-1.04800";
var numb: [0..1] int;
var n = arr.split(".");
numb[0] = abs(n[0]:int);
numb[1] = n[1]:int;
writeln((1.0 * numb[1]) / 10**(numb[1]:string).size);
writeln(numb[0]);
var n1: real = ((1.0 * numb[1]) / 10**(numb[1]:string).size) + numb[0];
writeln(n1);
if n[0][0] == '-' then n1 = n1 * (-1.0);
writeln(n1);
