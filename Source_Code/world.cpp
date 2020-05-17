//  Adam Sokolow 
//  adam.sokolow@duke.edu 
//  Dust simulation for Dr. Sen 
//  Dated July 29 2005

/*	Edited by Kevin VanSlyke
kgvansly@buffalo.edu
Dated Jan 2 2016	*/

#include <ctime> 
#include <cstdlib> 
#include <cmath> 
#include "world.h" 
#include <malloc.h> 
#include <stdio.h>  // for using fopen 

//Default constructor
world::world()
{
	myXSites = 500;
	myYSites = 500;
	myXMom = 10;
	myYMom = 10;
	myNegYMom = 0;	
	
	myWorld = (int**)malloc(myYSites*sizeof(int*));
	for (int c = 0; c < myYSites; ++c)
	{
		myWorld[c] = (int *)malloc(myXSites*sizeof(int));
		for (int d = 0; d < myXSites; ++d)
			myWorld[c][d] = -1;
	}
	myList = new dust_list(myWorld);
	myList->setMaxXLoc(myXSites);
	myList->setMaxYLoc(myYSites);
	myList->setMaxXMom(myXMom);
	myList->setMaxYMom(myYMom);
	myList->setNegYMom(myNegYMom);
}

//Constructor with size and momentum inputs
world::world(int x, int y, int XMom, int YMom, int negYMom)
{
	myXSites = x;
	myYSites = y;

	myXMom = XMom;
	myYMom = YMom;
	myNegYMom = negYMom;
	myWorld = (int**)malloc(myYSites*sizeof(int*));
	for (int c = 0; c < myYSites; ++c)
	{
		myWorld[c] = (int *)malloc(myXSites*sizeof(int));
		for (int d = 0; d < myXSites; ++d)
			myWorld[c][d] = -1;
	}
	myList = new dust_list(myWorld);
	myList->setMaxXLoc(myXSites);
	myList->setMaxYLoc(myYSites);
	myList->setMaxXMom(myXMom);
	myList->setMaxYMom(myYMom);
	myList->setNegYMom(myNegYMom);
}

// Copy constructor 
world::world(const world  & w)
{
	for (int c = 0; c < myXSites; ++c)
	{
		free(myWorld[c]);
	}
	free(myWorld);
	delete myList;
	myXSites = myYSites = 0;

	myXSites = w.myXSites;
	myYSites = w.myYSites;

	myXMom = w.myXMom;
	myYMom = w.myYMom;
	myNegYMom = w.myNegYMom;

	myWorld = (int**)malloc(myYSites*sizeof(int*));
	for (int c = 0; c < myYSites; ++c)
	{
		myWorld[c] = (int *)malloc(myXSites*sizeof(int));
		for (int d = 0; d < myXSites; ++d)
			myWorld[c][d] = w.myWorld[c][d];
	}
	myList = w.myList;
}

// Destructor, frees memory 
world ::~world()
{
	for (int c = 0; c < myYSites; ++c)
	{
		free(myWorld[c]);
	}
	free(myWorld);

	myXMom = myYMom = myNegYMom = myXSites = myYSites = 0;
	delete myList;
}

const world  &
world ::operator = (const world  & rhs)
{
	if (this != &rhs)                           // don't assign to self! 
	{
		for (int c = 0; c < myXSites; ++c)
		{
			free(myWorld[c]);
		}
		free(myWorld);
		delete myList;

		myXMom = myYMom = myNegYMom = myXSites = myYSites = 0;

		myXSites = rhs.myXSites;
		myYSites = rhs.myYSites;

		myList = rhs.myList;

		myWorld = (int**)malloc(myYSites*sizeof(int*));
		for (int c = 0; c < myYSites; ++c)
		{
			myWorld[c] = (int *)malloc(myXSites*sizeof(int));
			for (int d = 0; d < myXSites; ++d)
				myWorld[c][d] = rhs.myWorld[c][d];
		}
	}
	return *this;
}

void world::setWorld(int x, int y, int id)
{
	myWorld[y][x] = id;
}

void world::setWorld(std::vector <int> X, std::vector <int> Y, int id)
{
	for (unsigned int c = 0; c < X.size(); c++)
	{
		myWorld[Y[c]][X[c]] = id;
	}
}

int world::grainNumAt(int x, int y)
{
	if ( x < myXSites && y < myYSites)
		return myWorld[y][x];
	else
	{
		return -1;
		std::cout << "Error occured: attempted fetching grain number at out of bound location." << std::endl;
	}
		
}

int world::getMaxXSize()
{
	return myXSites;
}

int world::getMaxYSize()
{
	return myYSites;
}

//Adds dust grains one at a time until the number of grians specified by a parameter is reached.
void world::populateWorld(int numDust, int low, int high)
{
	std::cout << "Populating World!" << std::endl;
	while (myList->getTotal() < numDust)
	{
//		std::cout << ".";
		myList->addGrain(low, high);
	}
	std::cout << std::endl;
}

//As above, but also adds a single filter
void world::populateWorld(int numDust, int low, int high, int filterGap, int filterWidth, int filterLength)
{
	std::cout << "Populating World!" << std::endl;

	myList->addGrain(filterGap, filterWidth, filterLength); //add first filter

	while (myList->getTotal() <= numDust)
	{
		myList->addGrain(low, high);
//		std::cout << "-";
	}
	std::cout << std::endl;

}

void world::bimodalPopulateWorld(int numDust, int low1, int high1, int low2, int high2, int filterGap, int filterWidth, int filterLength)
{
	std::cout << "Populating World!" << std::endl;

	myList->addGrain(filterGap, filterWidth, filterLength); //add first filter

	while (myList->getTotal() <= numDust/2)
	{
		myList->addGrain(low1, high1);
//		std::cout << "-";
	}
	while (myList->getTotal() <= numDust)
	{
		myList->addGrain(low2, high2);
//		std::cout << "-";
	}
	std::cout << std::endl;

}

//As above, but also adds a second filter in addition to the first
//Probably not needed, could just call the previous function again with different inputs
void world::populateWorld(int numDust, int low, int high, int filterGap, int filterWidth, int filterLength, int  filter2Gap, int filter2Width, int filter2Length)
{
	std::cout << "Populating World!" << std::endl;

	myList->addGrain(filterGap, filterWidth, filterLength); //add first filter
	myList->addGrain2(filter2Gap, filter2Width, filter2Length); //add second filter

	while (myList->getTotal() <= numDust + 1)
	{
		myList->addGrain(low, high);
//		std::cout << ".";
	}
	std::cout << std::endl;
}

/*	//Routine to create oscillations in the filter to study dust response
//TODO: Test and get working for one filter
void world::shakefilter(int filterGap, int filterWidth, int filterLength, float timeF)
{
	updateWorld();
	myList->addGrain2(filterGap, filterWidth, filterLength);

}	*/

//Moves the simulation forward one time step.
void world::takeStep()
{
	//	std::cout<<"1"<<std::endl; 
	myList->moveStep(myWorld);

}

void world::updateWorld()
{
	dust_grain temp;

	for (int n = 0; n < myList->getTotal(); ++n)
	{
		temp = myList->getGrainByVecLoc(n);
		for (int d = 0; d < temp.getSize(); ++d)
			setWorld(temp.getXatc(d), temp.getYatc(d), temp.getID());
	}
}

void world::setProcOutputFolder(std::string dirName)
{
	procOutputFolder = dirName;
}

//The following function can be used to track dust movement via an output text file.  
void world::writingDust()
{
	dust_grain temp;

	std::string dustBeltFile = procOutputFolder + "/dustBeltX.txt";
	FILE * pFileBelt;
	// OUTPUT: number, size, y-position in a belt [240, 250]
	pFileBelt = fopen(dustBeltFile.c_str(), "a");

	std::string beltCountFile = procOutputFolder + "/beltCountX.txt";
	FILE * pBeltCount;
	pBeltCount = fopen(beltCountFile.c_str(), "a");

	std::vector<int> belt_dust;
	
	for (int n = 0; n < myList->getTotal(); ++n)
	{
		temp = myList->getGrainByVecLoc(n);
		if (temp.getYatc(0) > 240 && temp.getYatc(0) < 250) // The belt can be turned on/off here 
		{
			if (n > 0)
			{
				//Commented out so that we don't make a massive data file which doesn't really help analysis, this is more for debugging.
				fprintf(pFileBelt, "%d %d %d \n", temp.getID(), temp.getSize(), temp.getYatc(0)); //Y = up or down
				belt_dust.push_back(temp.getID());   // Add item at end of "belt" vector              
			}
		}
	}    // after done for "myTotal" (all grains) the "for" loop is closed 

	int belt_dustS;  // the size of this vector is the number of particles in the belt 
	belt_dustS = belt_dust.size();
	fprintf(pBeltCount, "%i \n", belt_dustS);
	belt_dust.clear();

	fclose(pFileBelt);
	fclose(pBeltCount);
}

//Outputs the number and location of overlapping particles.
void world::overlapingDust()
{
	//at t=ti     
	//check if grain1(i) = grain2(j)
	//write overlapinin in grains, 1,2     
	dust_grain tempi;
	dust_grain tempj;

	std::string dustOverlapFile = procOutputFolder + "/dustfileOverlap.txt";
	FILE * pFileOverlap;
	pFileOverlap = fopen(dustOverlapFile.c_str(), "a"); // overlaping

	for (int i = 0; i < myList->getTotal(); ++i)
	{
		for (int j = 0; j < myList->getTotal(); ++j)
		{
			tempi = myList->getGrainByVecLoc(i);
			tempj = myList->getGrainByVecLoc(j);
			if (i != j)
			{
				if ((tempi.getXatc(0) == tempj.getXatc(0)) && (tempi.getYatc(0) == tempj.getYatc(0)))
				{
					fprintf(pFileOverlap, "%d %d \n", i, j);
				}
			}
		}
	}
	fclose(pFileOverlap);
}

int ** world::getWorldArray()
{
	return myWorld;
}

int world::getCurNumDust()
{
	return myList->getTotal();
}
