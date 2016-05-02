// *******************************************************************
//  Random Generator for both Gaussian Distributions and Equal Probability
// *******************************************************************

//  Adam Sokolow
//  adam.sokolow@duke.edu
//  Dust simulation for Dr. Sen
//  Dated July 29 2005

#include <ctime>
#include <cstdlib>
#include <cmath>
#include "randomgen.h"

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>

//#include <mpi.h>



random_gen::random_gen()
{
	reseed();
}

random_gen::random_gen(int Seed)
{
	mySeed = Seed;
}

random_gen::random_gen(const random_gen  & gen)
{
	mySeed = gen.mySeed;
}


random_gen ::~random_gen ()
{

}

const random_gen  &
random_gen ::operator = (const random_gen  & rhs)
{
	if (this != &rhs)                           // don't assign to self!
	{
		mySeed = rhs.mySeed;
	}
	return *this;
}

int random_gen::getSeed()
{
	return mySeed;
}
double random_gen::Ran() // Ran() is uniform in 0..1
{
	const int IA = 16807, IC = 2147483647, IQ = 127773, IR = 2836;
	int h = mySeed / IQ;
	int l = mySeed % IQ;
	mySeed = IA * l - IR * h;

	if (mySeed <= 0)
	mySeed += IC;

	return mySeed / double(IC);
}
//***********gonsalves web
//Alejandro L. Garcia, Numerical Methods for Physics
//(Prentice Hall, Englewood Cliffs NJ, 2000)
double random_gen::rand() // rand() is uniform in 0..1
{
	const double a = 16807.0;
	const double m = 2147483647.0;
	double temp = a * mySeed;
	mySeed = (long)(fmod(temp,m));
	double rand = mySeed/m;
	return( rand );
}
//*********************ran1 Num Rec
double random_gen::ran1(long *idum) // ran1() is uniform in 0..1
{
	const int IA=16807, IM=2147483647;
	const int AM=(int)(1.0/IM), IQ=127773;
	const int IR=2836,NTAB=32;
	const int NDIV= (1+(IM-1)/NTAB);
	const double EPS=1.2e-7;
	const double RNMX=(1.0-EPS);
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if(*idum <= 0|| !iy)
	{
		if(-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for(j=NTAB+7; j>=0;j--)
		{
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum +=IM;
			if (j <NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum <0) *idum +=IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j]=*idum;
	if ((temp=(float)AM*iy) > RNMX) return RNMX;
	else return temp;
}
//*********************
double random_gen::Abs(double Nbr) //	return (Nbr >= 0) ? Nbr : -Nbr;
{
	if( Nbr >= 0 )
		return Nbr;
	else
		return -Nbr;
}
//This was also covered in PHY410
double random_gen::Gauss() // normal random variable generator
{
	double PI = 4.0*std::atan(1.0);
	double R1 = -1*std::log(1.0-rand());
	double R2 = 2.0*PI*rand();
	R1 = std::sqrt(2.0*R1);

	return (R1*std::cos(R2));
	//return R1*std::sin(R2);
}
/* boxmuller.c  //Gaussian distribution//////
Implements the Polar form of the Box-Muller
Transformation

c) Copyright 1994, Everett F. Carter Jr.
Permission is granted by the author to use
this software for any application provided this
copyright notice is preserved.

*/
float random_gen::Gauss2()	// normal random variable generator
{				            // mean m, standard deviation s
	float x1, x2, w, y1;
	static float y2;
	static int use_last = 0;

	if (use_last)		    // use value from previous call
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do 
		{
			x1 = (float)(2.0 * rand() - 1.0);
			x2 = (float)(2.0 * (float)rand() - 1.0);
			w = x1 * x1 + x2 * x2;
		} 
		while ( w >= 1.0 );

		w = (float)sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}
	return (float)Abs( 0.5 + y1 * 0.1 ); // NUMERIC: mean value & standard deviation
}

///////////////////////////////

//***********gonsalves web
//Alejandro L. Garcia, Numerical Methods for Physics
//(Prentice Hall, Englewood Cliffs NJ, 2000)
// Random number generator; Normal (Gaussian) dist.
double random_gen::randn()
{
	double randn = sqrt( -2.0*log(1.0 - rand()) ) * cos( 6.283185307 * rand() );
	return(randn);
}
//This can be found in Numerical Recipes
//It was also covered in PHY410 at SUNY Buffalo Computational Physics with Dr. Gonsalves.
void random_gen::reseed()
{
//    int npes, iam, ierr;
//   ierr = MPI_Comm_rank(MPI_COMM_WORLD, &iam);  /* get current process id */
//   ierr = MPI_Comm_size(MPI_COMM_WORLD, &npes);     /* get number of processes */
	std::time_t t;                   // usually an unsigned long int
	std::time(&t);                   // get seconds since Jan 1, 1900
	std::tm* pt = std::gmtime(&t);   // Universal (Greenwich Mean) Time
	mySeed = /*(iam/npes + iam%npes + iam) +*/ pt->tm_year + 70 * (pt->tm_mon + 12 * (pt->tm_hour +	23 * (pt->tm_min + 59 * pt->tm_sec)));
}