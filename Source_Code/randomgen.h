	
//		Adam	Sokolow	
//		adam.sokolow@duke.edu	
//		Dust	simulation	for	Dr.	Sen	
//		Dated	July	29	2005	
	
#ifndef	RANDOM_GEN_H	
	
	
#include<iostream>	
#include<cstdlib>	
	
//	*******************************************************************	
//	Create	a	Random	Generator	Class	that	has	a	variety	of	Options	
//	As	in	creating	Random	numbers	with	Gaussian	distribution,	and	etc.	
//	*******************************************************************	
	
		
class	random_gen	
{	
	public:			
	//	constructors/destructor	
	random_gen(	);											//	default	constructor	(size==0)	
	random_gen(int Seed);									//	initial	size	of	vector	is	size	
	random_gen(const random_gen	& gen);					//	copy	constructor	
	~random_gen();											//	destructor	
	//	assignment	
	const	random_gen	&	operator	=	(	const	random_gen	&	gen	);	
	//	accessors	
	int	getSeed();											//	Return	current	seed	value	
	double Ran();											//Return	a	random	real	number	between	0	and	1	
	double Gauss();

	double rand();
	double ran1(long	*idum);
	double randn();	
	double Abs(double	Nbr);	
	float Gauss2();		
	//	modifiers	
	void reseed();										//	Force	Reseeding									
	private:																																						
	int	mySeed;											//Current	Seed	
};	
#define	RANDOM_GEN_H	
#endif	
