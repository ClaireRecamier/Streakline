#include "streakline.h"
#include "nr_rand.c"
#include "utility.c"
#include <stdio.h>


const double M_sun = 4 * pi *pi;
const double mau = 0.00000000000668458; //meters to AU
const double secyr = 60*60*24*365.24; //seconds to years
const double aukpc = 0.00000000484814;//au to kpc


// Simulation parameters
double Mcli, Mclf, Rcl, dt;

// potential definitions
int par_perpotential[11] = {1, 3, 6, 6, 11, 4, 15, 13, 14, 19, 26};

int stream(double *x0, double *v0, double *xm1, double *xm2, double *xm3, double *xp1, double *xp2, double *xp3, double *vm1, double *vm2, double *vm3, double *vp1, double *vp2, double *vp3, double *par, double *offset, int potential, int integrator, int N, int M, double mcli, double mclf, double rcl, double dt_, FILE *fpt)
{
	int i,j, k=0, Napar, Ne, imin=0; //k is number of stream particles released
	double x[3], v[3], xs[3], vs[3], omega[3], om, sign=1., back=-1., r, rp, rm, vtot, vlead, vtrail, dM, Mcl, dR, dRRj, time=0.; //x is array of position of cluster, treated as a point mass?
	double *xc1, *xc2, *xc3, *Rj, *dvl, *dvt;
	long s1=560;

    double xlmc[3], vlmc[3] = {262000, 465000, 56000};
    if (potential==6){
        for(j=0;j<3;j++){
            xlmc[j] = par[12+j];
//             printf("%e ", xlmc[j]);
        }
//         printf("\n");
    }

	// number of output particles
	Ne=ceil((float)N/(float)M); //N number of desired timestep, particle is released every Mth timestep
	xc1 = (double*) malloc(N*sizeof(double)); //x coordinate of position of cluster, times number of timesteps? will store each x coordinate of cluster after each timestep
	xc2 = (double*) malloc(N*sizeof(double));
	xc3 = (double*) malloc(N*sizeof(double));
	Rj = (double*) malloc(Ne*sizeof(double));
	dvl = (double*) malloc(Ne*sizeof(double)); //velocity offset for leading tail, at each timestep
	dvt = (double*) malloc(Ne*sizeof(double)); //velocity offset of for trailing tail, at each timestep

	// Initial position and velocities
	t2t(x0,x); //copy array x0 to x
	t2t(v0,v);

	// Initial mass and mass loss rate
	Mcli = mcli; //initial mass of cluster
	Mclf = mclf; //final mass of cluster
	Mcl = Mcli; //current mass of cluster?
	dM = (Mcli-Mclf)/(double)N; //tott2,al mass lost, divided by number of timesteps

	// Cluster size
	Rcl = rcl;
/*
	// Position offset
	dR = offset[0];

	// Initialize velocity offsets, drawn from a Maxwell distribution
	double r1,r2,r3; //x y and z components of velocity offset?
	for(i=0;i<Ne;i++){
		// Leading tail velocity offsets
		r1=gasdev(&s1); //generate a standard normally distributed random variable?
		r2=gasdev(&s1);
		r3=gasdev(&s1);
		dvl[i]=sqrt(r1*r1 + r2*r2 + r3*r3)*offset[1]/3.;

		// Trailing tail dtvelocity offsets
		r1=gasdev(&s1);
		r2=gasdev(&s1);
		r3=gasdev(&s1);
		dvt[i]=sqrt(r1*r1 + r2*r2 + r3*r3)*offset[1]/3.;
	} */

	// Set up actual potential parameters;
	Napar = par_perpotential[potential]; //index 0 in array par_potential = 1 (parameters per potential?)
	double apar[Napar], apar_aux[11]; //array called apar of size 1, array called apar_aux of size 11
	initpar(potential, par, apar); //(0,input parameter, )

    if(potential==6){
        for(i=0;i<11;i++)
            apar_aux[i] = apar[i];
    }

	// Integrator switch
	void (*pt2dostep)(double*, double*, double*, int, double, double) = NULL; //a pointer to a function with no return value?


	if(integrator==0){
		// Leapfrog
		pt2dostep=&dostep; //set pointer to function to point to dostep function. Dont actually call function?
	}
	else if(integrator==1){
		// Runge-Kutta
		pt2dostep=&dostep_rk;
	}

	// Time step
	dt = dt_;

	///////////////////////////////////////
	// Backward integration (cluster only)

	if(integrator==0){ //if leapfrog
		//printf("cluster pos before first halfstep: %0.10f,%0.10f,%0.10f\n",aukpc * x[0],aukpc * x[1],aukpc * x[2]);
		//printf("cluster vel before first halfstep: %0.10f,%0.10f,%0.10f\n",aukpc * v[0],aukpc * v[1],aukpc * v[2]);
		dostep1(x,v,apar,potential,dt,back); //halfstep backwards in only velocity. X is an array of 3 - referring to cluster as one pointmass?
    //fprintf(fpt,"%0.30f,%0.30f,%0.30f,%0.30f\n",aukpc * x[0],aukpc * x[1],aukpc * v[0],aukpc * v[1]);
    //printf("cluster pos after first halfstep: %f,%f,%f\n",aukpc * x[0],aukpc * x[1],aukpc * x[2]);
		//printf("cluster vel after first halfstep: %0.15f,%0.10f,%f\n",aukpc * v[0],aukpc * v[1],aukpc * v[2]);
    imin=1;
        time = time + dt*back; //set current time to halfstep back

        if(potential==6){
            dostep1(xlmc,vlmc,apar_aux,4,dt,back);
            for(j=0;j<3;j++){
                apar[12+j] = xlmc[j];
//                 printf("%e ", xlmc[i]);
            }
//             printf("\n");
        }
	}
	for(i=imin;i<N;i++){ //N is number of desired timesteps. if leapfrog, imin=1, if RK, imin=0
		(*pt2dostep)(x,v,apar,potential,dt,back); //step backwards in pos and vel
		//printf("%0.30f,%0.30f,%0.30f,%0.30f\n",aukpc * x[0],aukpc * x[1],aukpc * v[0],aukpc * v[1]);
		//fprintf(fpt,"%0.30f,%0.30f,%0.30f,%0.30f\n",aukpc * x[0],aukpc * x[1],aukpc * v[0],aukpc * v[1]);
        time = time + dt*back;

        if(potential==6){
            (*pt2dostep)(xlmc,vlmc,apar_aux,4,dt,back);
            for(j=0;j<3;j++){
                apar[12+j] = xlmc[j];
//                 printf("%e ", xlmc[j]);
            }
//             printf("%d\n", i);
        }
	}
	if(integrator==0){ //if leapfrog, take another halfstep back in just velocity
		dostep1(x,v,apar,potential,dt,sign);//dostep1(x,v,apar,potential,dt,back); //
		//fprintf(fpt,"%0.30f,%0.30f,%0.30f,%0.30f\n",aukpc * x[0],aukpc * x[1],aukpc * v[0],aukpc * v[1]);

				if(potential==6){
            dostep1(xlmc,vlmc,apar_aux,4,dt,back);
            for(i=0;i<3;i++){
                apar[12+i] = xlmc[i];
//                 printf("%e ", xlmc[i]);
            }
//             printf("\n");
        }
    }

//     printf("%e", time);

	////////////////////////////////////////////
	// Forward integration (cluster and stream)

	// Initial step for the leapfrog integrator
	if (integrator==0){ //halfstep forward in only velocity if leapfrog
		dostep1(x,v,apar,potential,dt,sign);
		//for(j=0;j<3;j++) //loop through each coordinate
			//x[j]=x[j]-dt*v[j]; //full step backwards in position?
		//printf("cluster pos after first halfstep: %f,%f,%f\n",x[0],x[1],x[2]);
		//printf("cluster vel after first halfstep: %f,%f,%f\n",v[0],v[1],v[2]);
		//fprintf(fpt,"%f,%f,%f,%f\n",aukpc * x[0],aukpc * x[1],aukpc * v[0],aukpc * v[1]);
        if(potential==6){
            dostep1(xlmc,vlmc,apar_aux,4,dt,sign);
            for(j=0;j<3;j++)
                xlmc[j]=xlmc[j]-dt*vlmc[j];
            for(j=0;j<3;j++)
                apar[12+j] = xlmc[j];
        }

		//dostep(x,v,apar,potential,dt,sign); //step forward in both pos and vel
		imin=1; //indicates how many steps out of N have already been completed

        if(potential==6){
            dostep(xlmc,vlmc,apar_aux,4,dt,sign);
            for(i=0;j<3;j++)
                apar[12+j] = xlmc[j];
        }

		// Update output arrays
		t2n(x, xc1, xc2, xc3, 0); //fill array xc1, xc2, xc3: xc1[0] = x[0], xc2[0] = x[1], xc3[0]=x[2]. Update array xc1, xc2, and xc3, which will store all the positions of the cluster at each timestel?
	}

	// Subsequent steps
	for(i=imin;i<N;i++){ //for each timestep
		//Mcl-=dM; //decrease current mass

		(*pt2dostep)(x,v,apar,potential,dt,sign); //move cluster forward
		//fprintf(fpt,"%f,%f,%f,%f\n",aukpc * x[0],aukpc * x[1],aukpc * v[0],aukpc * v[1]);
        if(potential==6){
            (*pt2dostep)(xlmc,vlmc,apar_aux,4,dt,sign);
            for(j=0;j<3;j++)
                apar[12+j] = xlmc[j];
        }

        // Store cluster position
		t2n(x, xc1, xc2, xc3, i); //store cluster position from x in xc1 xc2 xc3

		// Propagate previously released stream particles
		for(j=0;j<k;j++){ //for each previously released stream particle
			// Inner particle
			n2t(xs, xm1, xm2, xm3, j); //fill array xs with index j of xm1, xm2, and xm3
			n2t(vs, vm1, vm2, vm3, j); //fill array vs with j of vm1, vm2, and vm3
			//printf("%f,%f,%f",xm1[j],xm2[j],xm3[j]);
			dostep_stream(x,xs,vs,apar,potential,Mcl,dt,sign); //xs and vs contain updated positions and velocities of stream.
// 			(*pt2dostep)(xs,vs,apar,potential,dt,sign);

			// Update
			t2n(xs, xm1, xm2, xm3, j); //copy xs (position of cluster) into xm1[j], xm2[j], and xm3[j]
			t2n(vs, vm1, vm2, vm3, j);

			// Outer particle
			n2t(xs, xp1, xp2, xp3, j); //fill array xs with values from xp123
			n2t(vs, vp1, vp2, vp3, j);
			dostep_stream(x,xs,vs,apar,potential,Mcl,dt,sign);//advance particle forward one timestep
// 			(*pt2dostep)(xs,vs,apar,potential,dt,sign);

			// Update
			t2n(xs, xp1, xp2, xp3, j); //copy position of particle into xp1, xp2, xp3. This is the same value of j as for the inner particle, whats the point of updating these arrays if they're not used?
			t2n(vs, vp1, vp2, vp3, j);

			//fprintf(fpt,",%f,%f,%f,%f,%f,%f,%f,%f",xm1[j],xm2[j],vm1[j],vm2[j],xp1[j],xp2[j],vp1[j],vp2[j]);
			//printf("stars %i at step %i %f,%f,%f,%f,%f,%f,%f,%f",j,i,aukpc * xm1[j],aukpc * xm2[j],aukpc * vm1[j],aukpc * vm2[j],aukpc * xp1[j],aukpc * xp2[j],aukpc * vp1[j],aukpc * vp2[j]);

		}

		/*
				if(i%M==0){

					// Release only at every Mth timestep
					// Jacobi tidal radius
					Rj[k]=jacobi(x, v, apar, potential, Mcl);
					r=len(x); //distance of cluster from origin
					rm=(r-Rj[k])/r; //inner particle
					rp=(r+Rj[k])/r; //outer particle

					// Angular velocity of cluster
					omega[0]=x[1]*v[2]-x[2]*v[1]; // r x v
					omega[1]=x[2]*v[0]-x[0]*v[2];
					omega[2]=x[0]*v[1]-x[1]*v[0];
					om=len(omega)/(r*r); // magnitude of: r x v / r^2
					vtot=len(v); //magnitude of velocity of cluster
					vlead=(vtot-om*Rj[k])/vtot; //velocity of cluster, minus velocity of particle, divided by velocity of cluster?
					vtrail=(vtot+om*Rj[k])/vtot;

					dvl[k]/=r;
					dvt[k]/=r;

					// Generate 2 new stream particles at the tidal radius
					dRRj = dR*Rj[k];
					// Inner particle (leading tail)
					xm1[k]=x[0]*rm + dRRj*gasdev(&s1);
					xm2[k]=x[1]*rm + dRRj*gasdev(&s1);
					xm3[k]=x[2]*rm + dRRj*gasdev(&s1);
					vm1[k]=v[0]*vlead - dvl[k]*x[0];
					vm2[k]=v[1]*vlead - dvl[k]*x[1];
					vm3[k]=v[2]*vlead - dvl[k]*x[2];

					// Outer particle (trailing tail)
					xp1[k]=x[0]*rp + dRRj*gasdev(&s1);
					xp2[k]=x[1]*rp + dRRj*gasdev(&s1);
					xp3[k]=x[2]*rp + dRRj*gasdev(&s1);
					vp1[k]=v[0]*vtrail + dvt[k]*x[0];
					vp2[k]=v[1]*vtrail + dvt[k]*x[1];
					vp3[k]=v[2]*vtrail + dvt[k]*x[2];

					k++;
				}

				time = time + dt*sign;
			}
		*/

		if (i % M == 0) {
				//EJECT particle at Mth timestep
				r=len(x); //distance of cluster from origin
				omega[0]=x[1]*v[2]-x[2]*v[1]; // r x v
				omega[1]=x[2]*v[0]-x[0]*v[2];
				omega[2]=x[0]*v[1]-x[1]*v[0];
		    //printf("omega of cluster: %f,%f,%f\n",omega[0],omega[1],omega[2]);
				om=len(omega)/(r*r); // magnitude of: r x v / r^2
				Rj[k]=jacobi(x, v, apar, potential, Mcl);
				rm=(r-Rj[k])/r; //inner particle
				rp=(r+Rj[k])/r; //outer particle
				//printf("tidal radius: %f,%f,%f\n",Rj[k] * x[0]/r,Rj[k] * x[1]/r,Rj[k] * x[2]/r);
				vtot=len(v); //magnitude of velocity of cluster
				vlead=(vtot-(om*Rj[k]))/vtot; //velocity of cluster, minus velocity of particle, divided by velocity of cluster?
				vtrail=(vtot+(om*Rj[k]))/vtot;
				//printf("v val %f\n",vtot-om*0.25);
				// Inner particle (leading tail)

				xm1[k]=x[0]*rm;
				xm2[k]=x[1]*rm;
				xm3[k]=x[2]*rm;
				vm1[k]=v[0]*vlead;
				vm2[k]=v[1]*vlead;
				vm3[k]=v[2]*vlead;
				//printf("cluster pos at ejection: %f,%f,%f\n",x[0],x[1],x[2]);
		  //printf("cluster vel at ejection: %f,%f,%f\n",v[0],v[1],v[2]);
			//printf("x leading particle at ejection: %f,%f,%f\n",xm1[k],xm2[k],xm3[k]);
		  //printf("v leading particle at ejection: %f,%f,%f\n",vm1[k],vm2[k],vm3[k]);

				// Outer particle (trailing tail)
				xp1[k]=x[0]*rp;
				xp2[k]=x[1]*rp;
				xp3[k]=x[2]*rp;
				vp1[k]=v[0]*vtrail;
				vp2[k]=v[1]*vtrail;
				vp3[k]=v[2]*vtrail;
				//printf("x trailing particle at ejection: %f,%f,%f\n",xp1[k],xp2[k],xp3[k]);
				//printf("v trailing particle at ejection: %f,%f,%f\n",vp1[k],vp2[k],vp3[k]);
				//fprintf(fpt,"%0.20f,%0.20f,%0.20f,%0.20f,%0.20f,%0.20f,%0.20f,%0.20f\n",aukpc * xm1[k],aukpc * xm2[k],aukpc * vm1[k],aukpc * vm2[k],aukpc * xp1[k],aukpc * xp2[k],aukpc * vp1[k],aukpc * vp2[k]);
				k++;

				//printf("%f,%f,%f,%f\n",xm1[0],xm2[0],xp1[0],xp2[0]);
				//printf("%f,%f,%f\n",vm1[0],vm2[0],vm3[0]);
	}

		//fprintf(fpt,"\n");
}


//     printf("%e\n", time);

    if (integrator==0){ //final halfstep back in velocity if leapfrog
		dostep1(x,v,apar,potential,dt,back);
		fprintf(fpt,"%f,%f,%f,%f\n",aukpc * x[0],aukpc * x[1],aukpc * v[0],aukpc *v[1]);

		for(j=0;j<k;j++) {
			fprintf(fpt,"%f,%f,%f,%f,%f,%f,%f,%f\n",aukpc * xm1[j],aukpc * xm2[j],aukpc * vm1[j],aukpc * vm2[j],aukpc * xp1[j],aukpc * xp2[j],aukpc * vp1[j],aukpc * vp2[j]);
		}




        if(potential==6){
            dostep1(xlmc,vlmc,apar_aux,4,dt,back);
            for(i=0;i<3;i++){
                apar[12+i] = xlmc[i];
//                 printf("%e ", xlmc[i]);
            }
//             printf("\n");
        }
    }

	// Free memory
	free(xc1);
	free(xc2);
	free(xc3);
	free(Rj);
	free(dvl);
	free(dvt);

	return 0;
}

int orbit(double *x0, double *v0, double *x1, double *x2, double *x3, double *v1, double *v2, double *v3, double *par, int potential, int integrator, int N, double dt_, double direction)
{
	int i, Napar, imin=0;
	double x[3], v[3];

	// Initial position and velocities
	t2t(x0, x); //copy array x0 to v
	t2t(v0, v);

    // Set up actual potential parameters;
    //Napar = par_perpotential[potential];
	//double apar[Napar];
	//initpar(potential, par, apar);



	// Integrator switch
	void (*pt2dostep)(double*, double*, double*, int, double, double) = NULL; //is this a function or a variable?

	if(integrator==0){
		// Leapfrog
		pt2dostep=&dostep;
	}
	else if(integrator==1){
		// Runge-Kutta
		pt2dostep=&dostep_rk;
	}

	// Time step
	dt = dt_;

	///////////////////////////////////////
	// Orbit integration

	if(integrator==0){
		dostep1(x,v,par,potential,dt,direction);

		// Record
		t2n(x, x1, x2, x3, 0); //record initial position at index 0
		t2n(v, v1, v2, v3, 0); //record vel at t=0.5dt at index 0
		//printf("%f\n",x1[0]);
		imin=1;
	}
	for(i=imin;i<N;i++){
		(*pt2dostep)(x,v,par,potential,dt,direction);

		// Record
		t2n(x, x1, x2, x3, i);
		t2n(v, v1, v2, v3, i);
		//printf("%f\n",x1[i]);
	}
	if(integrator==0){
		dostep1(x,v,par,potential,dt,-1.);

		// Record
		t2n(x, x1, x2, x3, N-1);
		t2n(v, v1, v2, v3, N-1); //record vel at t=
	}

	/*
	printf("x1");
	for (int a=0;a<N;a++){
		printf("%f",x1[a]);
	}
	*/


	return 0;
}

void dostep(double *x, double *v, double *par, int potential, double deltat, double sign)
{	// Evolve point particle from x0, v0, for a time deltat in a given potential
	// evolve forward for sign=1, backwards for sign=-1
	// return final positions and velocities in x, v

	int i, j, Nstep;
	double xt[3], vt[3], at[3], dts;

	dts=sign*dt;			// Time step with a sign
	Nstep=(int) (deltat/dt);	// Number of steps to evolve

	for(i=0;i<Nstep;i++){
		// Forward the particle using the leapfrog integrator
		for(j=0;j<3;j++)
			xt[j]=x[j]+dts*v[j];
		force(xt, at, par, potential);


        for(j=0;j<3;j++)
			vt[j]=v[j]+dts*at[j];

		// Update input vectors to current values
		for(j=0;j<3;j++){
			x[j]=xt[j];
			v[j]=vt[j];
		}
	}

}

void dostep1(double *x, double *v, double *par, int potential, double deltat, double sign)
{	// Make first step to set up the leapfrog integration

	double a[3], dts;

	dts=sign*dt;
	force(x, a, par, potential);
	//printf("dt %f calc from halfstep: %f,%f,%f\n",dts,0.5*dts*a[0],0.5*dts*a[1],0.5*dts*a[2]);
  //printf("x from halfstep %f,%f,%f\n",x[0],x[1],x[2]);

  //printf("a from halfstep %0.30f,%0.30f,%0.30f\n",a[0],a[1],a[2]);
	v[0]=v[0]+0.5*dts*a[0];
	v[1]=v[1]+0.5*dts*a[1];
	v[2]=v[2]+0.5*dts*a[2];
}

void dostep_rk(double *x, double *v, double *par, int potential, double deltat, double sign)
{	// Evolve point particle from x0, v0, for a time deltat in a given potential
	// evolve forward for sign=1, backwards for sign=-1
	// return final positions and velocities in x, v
	// prototype for Runge-Kutta integrator

	int i;
	double xt1[3], xt2[3], xt3[3], vt1[3], vt2[3], vt3[3], a[3], at1[3], at2[3], at3[3], dts, dt2;

	dts=sign*dt;			// Time step with a sign
	dt2=dts/2.;

	// Initial values
	force(x, a, par, potential);

	// First half-step
	for(i=0;i<3;i++){
		xt1[i]=x[i]+dt2*v[i];
		vt1[i]=v[i]+dt2*a[i];
	}
	force(xt1,at1,par,potential);

	// Second half-step
	for(i=0;i<3;i++){
		xt2[i]=x[i]+dt2*vt1[i];
		vt2[i]=v[i]+dt2*at1[i];
	}
	force(xt2,at2,par,potential);

	// Third step
	for(i=0;i<3;i++){
		xt3[i]=x[i]+dts*vt2[i];
		vt3[i]=v[i]+dts*at2[i];
	}
	force(xt3,at3,par,potential);

	// Final Runge-Kutta evaluation
	for(i=0;i<3;i++){
		x[i]+=dts/6.*(v[i]+2.*(vt1[i]+vt2[i])+vt3[i]);
		v[i]+=dts/6.*(a[i]+2.*(at1[i]+at2[i])+at3[i]);
	}

}

void dostep_stream(double *xc, double *x, double *v, double *par, int potential, double Mcl, double deltat, double sign)
{	// Same as dostep, except that stream particles also feel the Plummer potential from a cluster
	//x and v store position/vel of particle, xc stores position of cluster
	int i, j, Nstep;
	double xt[3], vt[3], at[3], xr[3], ar[3], dts;

	dts=sign*dt;			// Time step with a sign
	Nstep=(int) (deltat/dt);	// Number of steps to evolve

	for(i=0;i<Nstep;i++){
		// Forward the particle using the leapfrog integrator
		for(j=0;j<3;j++){
			xt[j]=x[j]+dts*v[j]; //update position of particle
			xr[j]=xc[j]-xt[j]; //xr store position of particle relative to cluster. what is the purpose, its not being used?
		}
		force(xt, at, par, potential); //acelerations, based on xt, get stored in at
//         if(potential==7) par[12]+=dts;
        force_plummer(xr,ar,Mcl); //accelerations, based on xr, get stored in ar
				//printf("xt: %f, %f, %f, a: %f,%f,%f,plummer:,%f,%f,%f\n",xt[0],xt[1],xt[2],at[0],at[1],at[2],ar[0],ar[1],ar[2]);

		for(j=0;j<3;j++)
			vt[j]=v[j]+dts*(at[j]+ar[j]); //update velocity with both acceleration vectors

		// Update input vectors to current values
		for(j=0;j<3;j++){
			x[j]=xt[j];
			v[j]=vt[j];
		}
	}
}

void force(double *x, double *a, double *par, int potential)
{
	int i;
	double r, aux, aux2;

	if(potential==0){
		// Point mass potential
		// par = [Mtot]
		r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
		for(i=0;i<3;i++)
			a[i]=-par[0]*x[i]/(r*r*r);
		//printf("acceleration: {%f,%f}\n", a[0],a[1]);

	}else if(potential==1){
		// Logarithmic potential, as defined by Koposov et al. (2010)
		// par = [Vc, q, q^2]
		r=x[0]*x[0] + x[1]*x[1] + x[2]*x[2]/(par[2]);
		aux=-par[0]/r;

		a[0]=aux*x[0];
		a[1]=aux*x[1];
		a[2]=aux*x[2]/(par[1]);

	}else if(potential==2){
		// Triaxial logarithmic halo potential from Law & Majewski (2010)
		// par = [Vc^2, c1, c2, c3, c4, rhalo^2]
		r=par[1]*x[0]*x[0] + par[2]*x[1]*x[1] + par[3]*x[0]*x[1] + par[4]*x[2]*x[2] + par[5];
		aux=-par[0]/r;

		a[0]=aux*(2*par[1]*x[0] + par[3]*x[1]);
		a[1]=aux*(2*par[2]*x[1] + par[3]*x[0]);
		a[2]=aux*(2*par[4]*x[2]);

	}else if(potential==3){
		// Triaxial NFW halo potential, parameters similar to Law & Majewski (2010)
		// par = [GM, c1, c2, c3, c4, rhalo]
		r=sqrt(par[1]*x[0]*x[0] + par[2]*x[1]*x[1] + par[3]*x[0]*x[1] + par[4]*x[2]*x[2]);
		aux=0.5 * par[0] / (r*r*r) * (1./(1.+par[5]/r)-log(1.+r/par[5]));

		a[0]=aux*(2*par[1]*x[0] + par[3]*x[1]);
		a[1]=aux*(2*par[2]*x[1] + par[3]*x[0]);
		a[2]=aux*(2*par[4]*x[2]);
		//printf("acceleration: {%0.10f,%0.10f}\n", a[0],a[1]);

	}else if(potential==4){
		// Composite Galactic potential featuring a disk, bulge, and flattened NFW halo (from Johnston/Law/Majewski/Helmi)
		// par = [GMb, ab, GMd, ad, bd^2, GM, c1, c2, c3, c4, rhalo]

		//Hernquist bulge
		r=sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
		aux=-par[0]/(r * (r+par[1]) * (r+par[1]));

		a[0]=aux*x[0];
		a[1]=aux*x[1];
		a[2]=aux*x[2];

		//Miyamoto-Nagai disk
		aux2=sqrt(x[2]*x[2] + par[4]);
		r=sqrt(x[0]*x[0] + x[1]*x[1] + (par[3] + aux2) * (par[3] + aux2));
		aux=-par[2]/(r*r*r);

		a[0]+=aux*x[0];
		a[1]+=aux*x[1];
		a[2]+=aux*x[2]*(par[3] + aux2)/aux2;

		//Triaxial NFW Halo
		r=sqrt(par[6]*x[0]*x[0] + par[7]*x[1]*x[1] + par[8]*x[0]*x[1] + par[9]*x[2]*x[2]);
		aux=0.5 * par[5]/(r*r*r) * (1./(1.+par[10]/r)-log(1.+r/par[10]));

		a[0]+=aux*(2*par[6]*x[0] + par[8]*x[1]);
		a[1]+=aux*(2*par[7]*x[1] + par[8]*x[0]);
		a[2]+=aux*(2*par[9]*x[2]);

	}else if(potential==5){
		// Spherical NFW potential
		// par = [GM, Rh]
		r=sqrt(x[0]*x[0] + x[1]*x[1]*par[2] + x[2]*x[2]*par[3]);
		aux=par[0]/(r*r*r) * (1./(1.+par[1]/r)-log(1.+r/par[1]));

		a[0]=aux*x[0];
		a[1]=aux*x[1]*par[2];
		a[2]=aux*x[2]*par[3];
	}else if(potential==6){
        // Galactic potential + LMC
        // par = [GMb, ab, GMd, ad, bd^2, GM, q^2, rhalo, GMlmc, Xlmc, Ylmc, Zlmc]

        //Hernquist bulge
        r=sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
        aux=-par[0]/(r * (r+par[1]) * (r+par[1]));

        a[0]=aux*x[0];
        a[1]=aux*x[1];
        a[2]=aux*x[2];

        //Miyamoto disk
        aux2=sqrt(x[2]*x[2] + par[4]);
        r=sqrt(x[0]*x[0] + x[1]*x[1] + (par[3] + aux2) * (par[3] + aux2));
        aux=-par[2]/(r*r*r);

        a[0]+=aux*x[0];
        a[1]+=aux*x[1];
        a[2]+=aux*x[2]*(par[3] + aux2)/aux2;

        //Triaxial NFW Halo
        r=sqrt(par[6]*x[0]*x[0] + par[7]*x[1]*x[1] + par[8]*x[0]*x[1] + par[9]*x[2]*x[2]);
        aux=0.5 * par[5]/(r*r*r) * (1./(1.+par[10]/r)-log(1.+r/par[10]));

        a[0]+=aux*(2*par[6]*x[0] + par[8]*x[1]);
        a[1]+=aux*(2*par[7]*x[1] + par[8]*x[0]);
        a[2]+=aux*(2*par[9]*x[2]);

        // Point mass
        // added softening ~8pc, assuming X_LMC~-0.8kpc
        r=sqrt((x[0]-par[12])*(x[0]-par[12]) + (x[1]-par[13])*(x[1]-par[13]) + (x[2]-par[14])*(x[2]-par[14])) - 0.01*par[12];
        aux = par[11]/(r*r*r);

        a[0]+=aux*x[0];
        a[1]+=aux*x[1];
        a[2]+=aux*x[2];

    }else if(potential==7){
        // Galactic potential + LMC on a string
        // par = [GMb, ab, GMd, ad, bd^2, GM, q^2, rhalo, Mlmc, t]

        //Hernquist bulge
        r=sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
        aux=-par[0]/(r * (r+par[1]) * (r+par[1]));

        a[0]=aux*x[0];
        a[1]=aux*x[1];
        a[2]=aux*x[2];

        //Miyamoto disk
        aux2=sqrt(x[2]*x[2] + par[4]);
        r=sqrt(x[0]*x[0] + x[1]*x[1] + (par[3] + aux2) * (par[3] + aux2));
        aux=-par[2]/(r*r*r);

        a[0]+=aux*x[0];
        a[1]+=aux*x[1];
        a[2]+=aux*x[2]*(par[3] + aux2)/aux2;

        //Triaxial NFW Halo
        r=sqrt(par[6]*x[0]*x[0] + par[7]*x[1]*x[1] + par[8]*x[0]*x[1] + par[9]*x[2]*x[2]);
        aux=0.5 * par[5]/(r*r*r) * (1./(1.+par[10]/r)-log(1.+r/par[10]));

        a[0]+=aux*(2*par[6]*x[0] + par[8]*x[1]);
        a[1]+=aux*(2*par[7]*x[1] + par[8]*x[0]);
        a[2]+=aux*(2*par[9]*x[2]);

        //LMC on a string
        double xlmc[3]={-2.509654716638902e19, -1.2653311505262738e21, -8.319850498177284e20};
        double vlmc[3]={-57, -226, 221};
        r=sqrt((x[0]-xlmc[0]-vlmc[0]*par[12])*(x[0]-xlmc[0]-vlmc[0]*par[12]) + (x[1]-xlmc[1]-vlmc[1]*par[12])*(x[1]-xlmc[1]-vlmc[1]*par[12]) + (x[2]-xlmc[2]-vlmc[2]*par[12])*(x[2]-xlmc[2]-vlmc[2]*par[12]));
        aux=-par[11]/(r*r*r);

        a[0]+=aux*x[0];
        a[1]+=aux*x[1];
        a[2]+=aux*x[2];
    }else if(potential==8){
		// Composite Galactic potential featuring a disk, bulge, flattened NFW halo (from Johnston/Law/Majewski/Helmi) and perturbations from dipole expansion
		// par = [GMb, ab, GMd, ad, bd^2, GM, c1, c2, c3, c4, rhalo, a10, a11, a12]

		//Hernquist bulge
		r=sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
		aux=-par[0]/(r * (r+par[1]) * (r+par[1]));

		a[0]=aux*x[0];
		a[1]=aux*x[1];
		a[2]=aux*x[2];

		//Miyamoto-Nagai disk
		aux2=sqrt(x[2]*x[2] + par[4]);
		r=sqrt(x[0]*x[0] + x[1]*x[1] + (par[3] + aux2) * (par[3] + aux2));
		aux=-par[2]/(r*r*r);

		a[0]+=aux*x[0];
		a[1]+=aux*x[1];
		a[2]+=aux*x[2]*(par[3] + aux2)/aux2;

		//Triaxial NFW Halo
		r=sqrt(par[6]*x[0]*x[0] + par[7]*x[1]*x[1] + par[8]*x[0]*x[1] + par[9]*x[2]*x[2]);
		aux=0.5 * par[5]/(r*r*r) * (1./(1.+par[10]/r)-log(1.+r/par[10]));

		a[0]+=aux*(2*par[6]*x[0] + par[8]*x[1]);
		a[1]+=aux*(2*par[7]*x[1] + par[8]*x[0]);
		a[2]+=aux*(2*par[9]*x[2]);

        // Dipole moment
        a[0]+=par[13];
        a[1]+=par[11];
        a[2]+=par[12];
    }else if(potential==9){
		// Composite Galactic potential featuring a disk, bulge, flattened NFW halo (from Johnston/Law/Majewski/Helmi) and perturbations from dipole and quadrupole moment
		// par = [GMb, ab, GMd, ad, bd^2, GM, c1, c2, c3, c4, rhalo, a10, a11, a12, a20, a21, a22, a23, a24]

		//Hernquist bulge
		r=sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
		aux=-par[0]/(r * (r+par[1]) * (r+par[1]));

		a[0]=aux*x[0];
		a[1]=aux*x[1];
		a[2]=aux*x[2];

		//Miyamoto-Nagai disk
		aux2=sqrt(x[2]*x[2] + par[4]);
		r=sqrt(x[0]*x[0] + x[1]*x[1] + (par[3] + aux2) * (par[3] + aux2));
		aux=-par[2]/(r*r*r);

		a[0]+=aux*x[0];
		a[1]+=aux*x[1];
		a[2]+=aux*x[2]*(par[3] + aux2)/aux2;

		//Triaxial NFW Halo
		r=sqrt(par[6]*x[0]*x[0] + par[7]*x[1]*x[1] + par[8]*x[0]*x[1] + par[9]*x[2]*x[2]);
		aux=0.5 * par[5]/(r*r*r) * (1./(1.+par[10]/r)-log(1.+r/par[10]));

		a[0]+=aux*(2*par[6]*x[0] + par[8]*x[1]);
		a[1]+=aux*(2*par[7]*x[1] + par[8]*x[0]);
		a[2]+=aux*(2*par[9]*x[2]);

        // Dipole moment
        a[0]+=par[13];
        a[1]+=par[11];
        a[2]+=par[12];

        // Quadrupole moment
        a[0]+= x[0]*(par[18] - par[16]) + x[1]*par[14] + x[2]*par[17];
        a[1]+= x[0]*par[14] - x[1]*(par[16] + par[18]) + x[2]*par[15];
        a[2]+= x[0]*par[17] + x[1]*par[15] + x[2]*2*par[16];
    }else if(potential==10){
		// Composite Galactic potential featuring a disk, bulge, flattened NFW halo (from Johnston/Law/Majewski/Helmi) and perturbations from dipole, quadrupole and octupole moments
		// par = [GMb, ab, GMd, ad, bd^2, GM, c1, c2, c3, c4, rhalo, a10, a11, a12, a20, a21, a22, a23, a24, a30, a31, a32, a33, a34, a35, a36]

		//Hernquist bulge
		r=sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
		aux=-par[0]/(r * (r+par[1]) * (r+par[1]));

		a[0]=aux*x[0];
		a[1]=aux*x[1];
		a[2]=aux*x[2];

		//Miyamoto-Nagai disk
		aux2=sqrt(x[2]*x[2] + par[4]);
		r=sqrt(x[0]*x[0] + x[1]*x[1] + (par[3] + aux2) * (par[3] + aux2));
		aux=-par[2]/(r*r*r);

		a[0]+=aux*x[0];
		a[1]+=aux*x[1];
		a[2]+=aux*x[2]*(par[3] + aux2)/aux2;

		//Triaxial NFW Halo
		r=sqrt(par[6]*x[0]*x[0] + par[7]*x[1]*x[1] + par[8]*x[0]*x[1] + par[9]*x[2]*x[2]);
		aux=0.5 * par[5]/(r*r*r) * (1./(1.+par[10]/r)-log(1.+r/par[10]));

		a[0]+=aux*(2*par[6]*x[0] + par[8]*x[1]);
		a[1]+=aux*(2*par[7]*x[1] + par[8]*x[0]);
		a[2]+=aux*(2*par[9]*x[2]);

        // Dipole moment
        a[0]+=par[13];
        a[1]+=par[11];
        a[2]+=par[12];

        // Quadrupole moment
        a[0]+= x[0]*(par[18] - par[16]) + x[1]*par[14] + x[2]*par[17];
        a[1]+= x[0]*par[14] - x[1]*(par[16] + par[18]) + x[2]*par[15];
        a[2]+= x[0]*par[17] + x[1]*par[15] + x[2]*2*par[16];

        // Octupole moment
        a[0]+= par[19]*6.*x[0]*x[1] + par[20]*x[1]*x[2] + par[21]*(-2.*x[0]*x[1]) + par[22]*(-6.*x[0]*x[2]) + par[23]*(4.*x[2]*x[2] - x[1]*x[1] - 3.*x[0]*x[0]) + par[24]*2.*x[0]*x[2] + par[25]*3.*(x[0]*x[0] - x[1]*x[1]);
        a[1]+= par[19]*3.*(x[0]*x[0] - x[1]*x[1]) + par[20]*x[0]*x[2] + par[21]*(4.*x[2]*x[2] - x[0]*x[0] - 3.*x[1]*x[1]) + par[22]*(-6.*x[1]*x[2]) + par[23]*(-2.*x[0]*x[1]) + par[24]*(-2.*x[1]*x[2]) + par[25]*(-6.*x[0]*x[1]);
        a[2]+= par[20]*x[0]*x[1] + par[21]*8.*x[1]*x[2] + par[22]*(6.*x[2]*x[2] - 3.*x[0]*x[0] - 3.*x[1]*x[1]) + par[23]*8*x[0]*x[2] + par[24]*(x[0]*x[0] - x[1]*x[1]);
    }
}

void force_plummer(double *x, double *a, double Mcl)
{	// Calculate acceleration a at a position x from a cluster with a Plummer profile
	// Assumes global definitions of cluster mass Mcl and radius Rcl
	int i;
	double r, raux;

	r=len(x); //magnitude of distance btwn cluster and particle
	raux=pow(r*r+Rcl*Rcl, 1.5); //r squared plus plummer radius squared

	for(i=0;i<3;i++)
		//a[i]=G*Mcl*x[i]/raux;
		a[i]=Mcl*x[i]/raux;
}

void initpar(int potential, double *par, double *apar)
{
	if(potential==1){
		// Logarithmic potential, par = [Vc, q_phi]
		// apar = [Vc^2, q, q^2]
		apar[0]=par[0]*par[0];
		apar[1]=par[1];
		apar[2]=par[1]*par[1];

	}else if(potential==2){
		// Triaxial halo potential from Law & Majewski (2010)
		// par = [Vc, phi, q_1, q_2, q_z, rhalo]
		// apar = [Vc^2, c1, c2, c3, c4, rhalo^2]
		double cosphi, sinphi;

		cosphi=cos(par[1]);
		sinphi=sin(par[1]);

		apar[0]=par[0]*par[0];
		apar[1]=cosphi*cosphi/(par[2]*par[2]) + sinphi*sinphi/(par[3]*par[3]);
		apar[2]=cosphi*cosphi/(par[3]*par[3]) + sinphi*sinphi/(par[2]*par[2]);
		apar[3]=2*sinphi*cosphi*(1/(par[2]*par[2]) - 1/(par[3]*par[3]));
		apar[4]=1/(par[4]*par[4]);
		apar[5]=par[5]*par[5];

	}else if(potential==3){
		// Triaxial NFW halo potential from Law & Majewski (2010)
		// par = [V, rhalo, phi, q_1, q_2, q_z]
		// apar = [GM, c1, c2, c3, c4, rhalo]
		double cosphi, sinphi;

		cosphi=cos(par[2]);
		sinphi=sin(par[2]);

// 		apar[0]=G*Msun*pow(10,par[0]);
		apar[0]=par[0]*par[0]*par[1];
		apar[1]=cosphi*cosphi/(par[3]*par[3]) + sinphi*sinphi/(par[4]*par[4]);
		apar[2]=cosphi*cosphi/(par[4]*par[4]) + sinphi*sinphi/(par[3]*par[3]);
		apar[3]=2*sinphi*cosphi*(1/(par[3]*par[3]) - 1/(par[4]*par[4]));
		apar[4]=1/(par[5]*par[5]);
		apar[5]=par[1];
		//printf("initialization of parameters: %f,%f,%f,%f,%f,%f\n",apar[0],apar[1],apar[2],apar[3],apar[4],apar[5]);

	}else if(potential==0){
		// Point mass potential, par = [Mtot]
		// apar = [G*Mtot]
		apar[0]= par[0];

	}else if(potential==4){
		// Composite Galactic potential featuring a disk, bulge, and triaxial NFW halo (from Johnston/Law/Majewski/Helmi)
		// par = [GMb, ab, GMd, ad, bd, V, rhalo, phi, q_1, q_2, q_z]
		// apar = [GMb, ab, GMd, ad, bd^2, GM, c1, c2, c3, c4, rhalo]
		double cosphi, sinphi; //, tq, tphi;

		apar[0]=G*par[0];
		apar[1]=par[1];
		apar[2]=G*par[2];
		apar[3]=par[3];
		apar[4]=par[4]*par[4];

		cosphi=cos(par[7]);
		sinphi=sin(par[7]);

		apar[5]=par[5]*par[5]*par[6];
		apar[6]=cosphi*cosphi/(par[8]*par[8]) + sinphi*sinphi/(par[9]*par[9]);
		apar[7]=cosphi*cosphi/(par[9]*par[9]) + sinphi*sinphi/(par[8]*par[8]);
		apar[8]=2*sinphi*cosphi*(1/(par[8]*par[8]) - 1/(par[9]*par[9]));
		apar[9]=1/(par[10]*par[10]);
		apar[10]=par[6];

    }else if(potential==5){
		apar[0]=G*Msun*pow(10,par[0]);
		apar[1]=par[1];
		apar[2]=par[5]*par[5];
		apar[3]=par[6]*par[6];

	}else if(potential==6){
        // Galactic potential + LMC
        // par = [GMb, ab, GMd, ad, bd, Vh, rhalo, phi, q1, q2, qz, Mlmc, Xlmc, Ylmc, Zlmc]
        // apar = [GMb, ab, GMd, ad, bd^2, GM, c1, c2, c3, c4, rhalo, GMlmc, Xlmc, Ylmc, Zlmc]
        double cosphi, sinphi;

        apar[0]=G*par[0];
        apar[1]=par[1];

        apar[2]=G*par[2];
        apar[3]=par[3];
        apar[4]=par[4]*par[4];

        cosphi=cos(par[7]);
        sinphi=sin(par[7]);
        apar[5]=par[5]*par[5]*par[6];
        apar[6]=cosphi*cosphi/(par[8]*par[8]) + sinphi*sinphi/(par[9]*par[9]);
        apar[7]=cosphi*cosphi/(par[9]*par[9]) + sinphi*sinphi/(par[8]*par[8]);
        apar[8]=2*sinphi*cosphi*(1/(par[8]*par[8]) - 1/(par[9]*par[9]));
        apar[9]=1/(par[10]*par[10]);
        apar[10]=par[6];

        apar[11]=G*par[11];
        apar[12]=par[12];
        apar[13]=par[13];
        apar[14]=par[14];
    }else if(potential==7){
        // Galactic potential + LMC on a string
        // par = [GMb, ab, GMd, ad, bd, V, rhalo, phi, q_1, q_2, q_z, Mlmc]
        // apar = [GMb, ab, GMd, ad, bd^2, GM, c1, c2, c3, c4, rhalo, Mlmc, t]
        double cosphi, sinphi; //, tq, tphi;

        apar[0]=G*par[0];
        apar[1]=par[1];
        apar[2]=G*par[2];
        apar[3]=par[3];
        apar[4]=par[4]*par[4];

        cosphi=cos(par[7]);
        sinphi=sin(par[7]);
        apar[5]=par[5]*par[5]*par[6];
        apar[6]=cosphi*cosphi/(par[8]*par[8]) + sinphi*sinphi/(par[9]*par[9]);
        apar[7]=cosphi*cosphi/(par[9]*par[9]) + sinphi*sinphi/(par[8]*par[8]);
        apar[8]=2*sinphi*cosphi*(1/(par[8]*par[8]) - 1/(par[9]*par[9]));
        apar[9]=1/(par[10]*par[10]);
        apar[10]=par[6];
        apar[11]=G*par[11];
//         printf("%e\n", par[11]);
        apar[12]=0.;
    }else if(potential==8){
		// Composite Galactic potential featuring a disk, bulge, and triaxial NFW halo (from Johnston/Law/Majewski/Helmi)
		// par = [GMb, ab, GMd, ad, bd, V, rhalo, phi, q_1, q_2, q_z, a10, a11, a12]
		// apar = [GMb, ab, GMd, ad, bd^2, GM, c1, c2, c3, c4, rhalo, fa10, fa11, fa12]
		double cosphi, sinphi, f; //, tq, tphi;

		apar[0]=G*par[0];
		apar[1]=par[1];
		apar[2]=G*par[2];
		apar[3]=par[3];
		apar[4]=par[4]*par[4];

		cosphi=cos(par[7]);
		sinphi=sin(par[7]);

		apar[5]=par[5]*par[5]*par[6];
		apar[6]=cosphi*cosphi/(par[8]*par[8]) + sinphi*sinphi/(par[9]*par[9]);
		apar[7]=cosphi*cosphi/(par[9]*par[9]) + sinphi*sinphi/(par[8]*par[8]);
		apar[8]=2*sinphi*cosphi*(1/(par[8]*par[8]) - 1/(par[9]*par[9]));
		apar[9]=1/(par[10]*par[10]);
		apar[10]=par[6];

        f = sqrt(3./(4*pi));
        apar[11] = f*par[11];
        apar[12] = f*par[12];
        apar[13] = f*par[13];
    }else if(potential==9){
		// Composite Galactic potential featuring a disk, bulge, and triaxial NFW halo (from Johnston/Law/Majewski/Helmi)
		// par = [GMb, ab, GMd, ad, bd, V, rhalo, phi, q_1, q_2, q_z, a10, a11, a12, a20, a21, a22, a23, a24]
		// apar = [GMb, ab, GMd, ad, bd^2, GM, c1, c2, c3, c4, rhalo, fa10, fa11, fa12, fa20, fa21, fa22, fa23, fa24]
		double cosphi, sinphi, f; //, tq, tphi;

		apar[0]=G*par[0];
		apar[1]=par[1];
		apar[2]=G*par[2];
		apar[3]=par[3];
		apar[4]=par[4]*par[4];

		cosphi=cos(par[7]);
		sinphi=sin(par[7]);

		apar[5]=par[5]*par[5]*par[6];
		apar[6]=cosphi*cosphi/(par[8]*par[8]) + sinphi*sinphi/(par[9]*par[9]);
		apar[7]=cosphi*cosphi/(par[9]*par[9]) + sinphi*sinphi/(par[8]*par[8]);
		apar[8]=2*sinphi*cosphi*(1/(par[8]*par[8]) - 1/(par[9]*par[9]));
		apar[9]=1/(par[10]*par[10]);
		apar[10]=par[6];

        f = sqrt(3./(4*pi));
        apar[11] = f*par[11];
        apar[12] = f*par[12];
        apar[13] = f*par[13];

        f = 0.5*sqrt(15./pi);
        apar[14] = f*par[14];
        apar[15] = f*par[15];
        apar[16] = f/sqrt(3.)*par[16];
        apar[17] = f*par[17];
        apar[18] = f*par[18];
    }else if(potential==10){
        // Composite Galactic potential featuring a disk, bulge, and triaxial NFW halo (from Johnston/Law/Majewski/Helmi)
		// par = [GMb, ab, GMd, ad, bd, V, rhalo, phi, q_1, q_2, q_z, a10, a11, a12, a20, a21, a22, a23, a24, a30, a31, a32, a33, a34, a35, a36]
		// apar = [GMb, ab, GMd, ad, bd^2, GM, c1, c2, c3, c4, rhalo, fa10, fa11, fa12, fa20, fa21, fa22, fa23, fa24, fa30, fa31, fa32, fa33, fa34, fa35, fa36]
		double cosphi, sinphi, f; //, tq, tphi;

		apar[0]=G*par[0];
		apar[1]=par[1];
		apar[2]=G*par[2];
		apar[3]=par[3];
		apar[4]=par[4]*par[4];

		cosphi=cos(par[7]);
		sinphi=sin(par[7]);

		apar[5]=par[5]*par[5]*par[6];
		apar[6]=cosphi*cosphi/(par[8]*par[8]) + sinphi*sinphi/(par[9]*par[9]);
		apar[7]=cosphi*cosphi/(par[9]*par[9]) + sinphi*sinphi/(par[8]*par[8]);
		apar[8]=2*sinphi*cosphi*(1/(par[8]*par[8]) - 1/(par[9]*par[9]));
		apar[9]=1/(par[10]*par[10]);
		apar[10]=par[6];

        // dipole
        f = sqrt(3./(4*pi));
        apar[11] = f*par[11];
        apar[12] = f*par[12];
        apar[13] = f*par[13];

        // quadrupole
        f = 0.5*sqrt(15./pi);
        apar[14] = f*par[14];
        apar[15] = f*par[15];
        apar[16] = f/sqrt(3.)*par[16];
        apar[17] = f*par[17];
        apar[18] = f*par[18];

        // octupole
        f = 0.25*sqrt(35./(2.*pi));
        apar[19] = f*par[19];
        apar[25] = f*par[25];

        f = 0.25*sqrt(105./pi);
        apar[20] = f*par[20]*2.;
        apar[24] = f*par[24];

        f = 0.25*sqrt(21./(2.*pi));
        apar[21] = f*par[21];
        apar[23] = f*par[23];

        apar[22] = 0.25*sqrt(7./pi)*par[22];
    }
}

double jacobi(double *x, double *v, double *par, int potential, double Mcl)
{	// Jacobi radius of a cluster, aka tidal radius
	// at the position x, velocity v, and in a given potential

	int i;
	double R, om, dpot, delta, r;
	double omega[3], x1[3], x2[3], a1[3], a2[3], dx[3]; //, da[3];

	// Radial distance of cluster from center of galaxy
	r=len(x);

	// Angular velocity
	omega[0]=x[1]*v[2]-x[2]*v[1];
	omega[1]=x[2]*v[0]-x[0]*v[2];
	omega[2]=x[0]*v[1]-x[1]*v[0];
	om=len(omega)/(r*r);
	//printf("pos cl: %f,%f,%f, vel cl: %f,%f,%f, omega: %f\n",x[0],x[1],x[2],v[0],v[1],v[2],om);

	// Potential derivative
	delta = 0.02 * kpcau;//delta=0.02*kpc;
	for(i=0;i<3;i++){
		x1[i]=x[i]/r*(r-delta); //unit vector * r - delta./ position of cluster ,slightly back in time
		x2[i]=x[i]/r*(r+delta); //position of cluster, slightly forward in time
	}
	//printf("x1: %f,%f,%f\n",x1[0],x1[1],x1[2]);
	//printf("x2: %f,%f,%f\n",x2[0],x2[1],x2[2]);

	force(x1, a1, par, potential);
	force(x2, a2, par, potential);
	//printf("a1: %0.20f,%0.20f,%0.20f\n",a1[0],a1[1],a1[2]);
	//printf("a2: %0.20f,%0.20f,%0.20f\n",a2[0],a2[1],a2[2]);
	//printf("len a1 %0.30f\n",len(a1));
	//printf("len a2 %0.30f\n",len(a2));
	for(i=0;i<3;i++){
		dx[i]=x1[i]-x2[i];
	}
	dpot=(len(a1)-len(a2))/len(dx);
	//printf("len x1-x2 %0.30f\n", len(dx));
	//printf("len a1-a2 %0.30f\n",len(a1)-len(a2));
	//printf("dpot: %0.30f,%0.30f,%0.30f\n",dpot*x[0]/r,dpot*x[1]/r,dpot*x[2]/r);

	// Jacobi radius
	//printf("mcl %0.10f\n",Mcl);
	R=pow(Mcl/fabs(om*om+dpot),1./3.);//R=pow(G*Mcl/fabs(om*om+dpot),1./3.); //fabs = absolute value of floating point number

	return R;
}


int main (void) {

	int N=6000;//6000
	double M = 1;
	int potential = 0;
	int integrator = 0;
  Rcl = 20 * 0.001 * kpcau;
	double x0[3] = {50.0*kpcau,0,0}; //initial positions
	double v0[3] = {0,0.5 * sqrt(1000000000 * M_sun/len(x0)),0}; //initial vel
	//double par[6] = {430.0 * pow(10,3) * mau * secyr, 19.5 * kpcau, 88.0, 0.855, 1.0, 1.2};//double par[1] = {M_sun};
	double par[1] = {1000000000 * M_sun};
	double dt = 1000000;
	int sign = 1;
	double *offset;
	double Ne = N/M;
	double *x1 = (double *)malloc(sizeof(double) * N);
	double *x2 = (double *)malloc(sizeof(double) * N);
	double *x3 = (double *)malloc(sizeof(double) * N);
	double *v1 = (double *)malloc(sizeof(double) * N);
	double *v2 = (double *)malloc(sizeof(double) * N);
	double *v3 = (double *)malloc(sizeof(double) * N);
	double *xm1 = (double *)malloc(sizeof(double) * Ne); //position of leading trail
	double *xm2 = (double *)malloc(sizeof(double) * Ne);
	double *xm3 = (double *)malloc(sizeof(double) * Ne);
	double *vm1 = (double *)malloc(sizeof(double) * Ne);
	double *vm2 = (double *)malloc(sizeof(double) * Ne);
	double *vm3 = (double *)malloc(sizeof(double) * Ne);
	double *xp1 = (double *)malloc(sizeof(double) * Ne); //position of trailing trail
	double *xp2 = (double *)malloc(sizeof(double) * Ne);
	double *xp3 = (double *)malloc(sizeof(double) * Ne);
	double *vp1 = (double *)malloc(sizeof(double) * Ne);
	double *vp2 = (double *)malloc(sizeof(double) * Ne);
	double *vp3 = (double *)malloc(sizeof(double) * Ne);

	FILE *fpt = fopen("ctest1.csv", "w+");
	fprintf(fpt,"xcl,ycl,xvcl,yvcl,xlt,ylt,xvlt,yvlt,xtt,ytt,xvtt,yvtt\n");

	stream(x0, v0, xm1, xm2,xm3, xp1, xp2, xp3, vm1, vm2, vm3, vp1, vp2, vp3, par, offset, potential, 0, N, M, 20000*M_sun,20000*M_sun, Rcl, dt, fpt);
	//orbit(x0, v0, x1, x2, x3, v1, v2, v3, par, potential, integrator, N, dt, sign);

	fclose(fpt);


}
