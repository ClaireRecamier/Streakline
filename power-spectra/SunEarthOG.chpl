//define constants pi, solarMass, number of timesteps
use IO;
config const iterations = 100, timestep = 0.01; //years
const pi = 3.141592653589793, daysPerYear = 365.24, solarMass = 4 * pi**2;

//define class/record of a body
record body {
  var position: 3*real;
  var velocity: 3*real;
  var mass: real;
}
var chfile = open("SEog.csv",iomode.cw); //create test.csv and open
var WritingChannel = chfile.writer(); //open writing channel to test.csv


//define an array of bodies
var bodies = [
  //SUN
  new body(position = (0.0,0.0,0.0),
          velocity = (0.0,0.0,0.0),
          mass = solarMass),

  /* earth */
  new body(position = ( 1.0,0.0,0.0),
           velocity = (0.0,sqrt(solarMass),0.0),
          mass =   (1.0/333000.0) * solarMass)
];
const numBodies = bodies.size;

//procedure to initialize velocity of sun. Checked, function works
proc initSun() {
  var p: 3*real = (0.0,0.0,0.0);
  for body in bodies do { //shortcut: could use reduce instead to shorten this function
    p += body.mass * body.velocity;
  }
  //const p = + reduce (for b in bodies do (b.velocity * b.mass));
  bodies[0].velocity = -p / solarMass;
  writeln("vel of sun:", bodies[0].velocity);
}

//procedure to advance by a timestep
proc advance(dt: real){
  //define array of accelerations
  var a: [0..numBodies-1] 3*real;
  force(a); //compute forces and store in accelerations array
  //WritingChannel.write(acc[1][0],",",acc[1][1],",",acc[1][2],",");

  for i in 0..numBodies-1 do{
    //add acceleration vector to velocity vector for new velocity
    bodies[i].velocity += dt * a[i];
    //add velocity vector to position for new position
    bodies[i].position += dt * bodies[i].velocity;
  }
  writeln(a[1]);
  writeln(bodies[1].velocity);
  WritingChannel.write(bodies[1].position[0]*bodies[1].velocity[0] + bodies[1].position[1]*bodies[1].velocity[1] + bodies[1].position[2]*bodies[1].velocity[2],",");
  WritingChannel.write(a[1][0]*bodies[1].velocity[0] + a[1][1]*bodies[1].velocity[1] + a[1][2]*bodies[1].velocity[2],",");
  WritingChannel.write(dist(a[1]),",");
  WritingChannel.write(dist(bodies[1].velocity),",");
  WritingChannel.write(dist(bodies[1].position),",");
}

//procedure to calculate energy
proc energy(){
  var energy: real = 0.0;
  //initialize relative position vector and its length as new variables
  var r: 3*real = (0.0,0.0,0.0);
  var r_len: real = 0.0;

  for i in 0..numBodies-1 do {
    //add KE
    energy += 0.5 * bodies[i].mass * (dist(bodies[i].velocity)**2);
    //subtract PE of interaction of this body with every other body (of which PE between has not already been calculated)
    for j in i+1..numBodies-1 do {
      //calculate relative position vector
      r = bodies[i].position - bodies[j].position;
      //calculate length of relative position vector
      r_len = dist(r);

      energy -= (bodies[i].mass * bodies[j].mass) / r_len;
    }
  }
  return energy;


}

//shortcut procedure to help calculate distance. Checked, function works.
proc dist((x,y,z)) {
  return sqrt(x**2 + y**2 + z**2);
}

proc force(acc){
  //define force vector
  var F: 3*real = (0.0,0.0,0.0);
  //initialize relative position vector and its length as new variables
  var r: 3*real;
  var r_len: real = 0.0;

  for i in 0..numBodies-1 do {
    for j in i+1..numBodies-1 do {
      //calculate relative position vector
      r = bodies[i].position - bodies[j].position;
      //calculate length of relative position vector
      r_len = dist(r);
      //Force btwn two bodies: F = G m1 m2 / r ^2 * unit vector
      F = (r/(r_len**3));
      //acc[i] -= F * bodies[j].mass;
      acc[j] += F * bodies[i].mass;
    }
  }
}

//write main procedure. Initialize velocity of sun, print energy. Advance by a timestep, print energy.
proc main() {
  //initSun();
  //writef("%.9r\n",energy());
  WritingChannel.write(bodies[1].position[0]*bodies[1].velocity[0] + bodies[1].position[1]*bodies[1].velocity[1] + bodies[1].position[2]*bodies[1].velocity[2],",");
  WritingChannel.write("0.0,0.0,");
  WritingChannel.write(dist(bodies[1].velocity),",");
  WritingChannel.write(dist(bodies[1].position),",");
  WritingChannel.write(bodies[1].position[0],",",bodies[1].position[1],",",bodies[1].position[2],",");
  WritingChannel.write(bodies[1].velocity[0],",",bodies[1].velocity[1],",",bodies[1].velocity[2],"\n");
  for i in 1..iterations do {
    advance(timestep);
    WritingChannel.write(bodies[1].position[0],",",bodies[1].position[1],",",bodies[1].position[2],",");
    WritingChannel.write(bodies[1].velocity[0],",",bodies[1].velocity[1],",",bodies[1].velocity[2],"\n");
    /*
    if i % 10000 == 0 { //print energy every 100th iteration
      writef("%.9r\n",energy());
    }
    */
  }
  WritingChannel.close();
  chfile.close();

  //writef("%.9r\n",energy());

  /* ALTERNATIVE main procedure, where change in energy gets printed instead of energy
  initSun();

  var init_energy = energy();
  var curr_energy = init_energy;
  var delta_energy = 0.0;
  writef("%.9r\n",init_energy);

  for i in 1..iterations do {
    advance(timestep);
    if i % 1000 == 0 { //print energy every 100th iteration
      curr_energy = energy();
      delta_energy = (curr_energy - init_energy)/init_energy;
      writef("%.9r\n",abs(delta_energy));
    }
  }
  //writef("%.9r\n",energy());
  */

}
