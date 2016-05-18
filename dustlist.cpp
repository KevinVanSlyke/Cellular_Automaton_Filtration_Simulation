//  Adam Sokolow
//  adam.sokolow@duke.edu
//  Dust simulation for Dr. Sen
//  Dated July 29 2005

/*	Edited by Kevin VanSlyke
kgvansly@buffalo.edu
Dated Jan 2 2016 */
// *******************************************************************
//  A List of Dust grains
// *******************************************************************

#include <ctime>
#include <cstdlib>
#include "dustgrain.h"
#include "randomgen.h"
#include "dustlist.h"
#include <malloc.h>
#include <cmath>
#include <algorithm>

#include <iostream>
#include <stdio.h>

//Default constructor
dust_list::dust_list()
{
	maxXLoc = 500;
	maxYLoc = 500;
	ptclsHandled = 0;
	enableSticking = true;
	enableMerging = false;
	enableSplitting = false;
	myTotal = 0;
	myGenerator = new random_gen();
	myDustList.resize(0);
	numTimeSteps = 0;
	uniqueID = 0;
}

dust_list::dust_list(int** world)
{
	maxXLoc = 500;
	maxYLoc = 500;
	ptclsHandled = 0;
	enableSticking = true;
	enableMerging = false;
	enableSplitting = false;
	myTotal = 0;
	myGenerator = new random_gen();
	myDustList.resize(0);
	refWorld = world;
	numTimeSteps = 0;
	uniqueID = 0;
}

//Constructor with parameters for particle creation.
dust_list::dust_list(int** world, int total, int low, int high)
{
	maxXLoc = 500;
	maxYLoc = 500;
	ptclsHandled = 0;
	enableSticking = true;
	enableMerging = false;
	enableSplitting = false;
	myTotal = 0;
	myGenerator = new random_gen();
	uniqueID = 0;
	myDustList.resize(total);
	for (int c = 0; c < total; ++c)
	{
		addGrain(low, high, newUniqueID());   //addGrain must create a dust of random size in a random open location
	}
	refWorld = world;
	numTimeSteps = 0;
}

//Copier
dust_list::dust_list(const dust_list  & dl)
{
	myGenerator = new random_gen();
	myTotal = dl.myTotal;
	myDustList.resize(myTotal);
	myDustList = dl.myDustList;
	refWorld = dl.refWorld;
	numTimeSteps = dl.numTimeSteps;
	uniqueID = dl.uniqueID;
	pillBoxes = dl.pillBoxes;
	dustDist = dl.dustDist;
	maxXLoc = dl.maxXLoc;
	maxYLoc = dl.maxYLoc;
	pBCounts = dl.pBCounts;
	poreJamTimer = dl.poreJamTimer;
	poreBlocked = dl.poreBlocked;
	potentialBlock = dl.potentialBlock;
	ptclsHandled = dl.ptclsHandled;
}

//Destructor
dust_list ::~dust_list()
{
	myTotal = 0;
	delete myGenerator;
	myGenerator = NULL;
	uniqueID = 0;
}

const dust_list  &
dust_list ::operator = (const dust_list  & rhs)
{
	if (this != &rhs)   // don't assign to self!
	{
		myGenerator = new random_gen();
		myTotal = rhs.myTotal;
		myDustList.resize(myTotal);
		myDustList = rhs.myDustList;
		numTimeSteps = rhs.numTimeSteps;
		refWorld = rhs.refWorld;
		uniqueID = rhs.uniqueID;
		pillBoxes = rhs.pillBoxes;
		dustDist = rhs.dustDist;
		maxXLoc = rhs.maxXLoc;
		maxYLoc = rhs.maxYLoc;
		pBCounts = rhs.pBCounts;
		poreJamTimer = rhs.poreJamTimer;
		poreBlocked = rhs.poreBlocked;
		potentialBlock = rhs.potentialBlock;
		ptclsHandled = rhs.ptclsHandled;

	}
	return *this;
}

int dust_list::getTotal()
{

	return myTotal;

}

//With conglomeration off this method checks that the new position is empty so that the current dust particle can move into the space.
bool dust_list::isOpen(int x, int y, int grain_id)
{
	if ((refWorld[y][x] == -1) || (refWorld[y][x] == grain_id))
		return true;
	else
		return false;
}

//Same as above but doesn't check if the new location overlaps with the current one.
bool dust_list::isOpen(int x, int y)
{
	if (refWorld[y][x] == -1)
		return true;
	else
		return false;
}

//Determines if cardinal neighboring elements are open for the current particle of size psz to be physically extended into
bool dust_list::isOpenSelf(std::vector<int> x, std::vector<int> y, int ranSite, int ranCardinal, int psz)
{
	int i, x1, y1;
	bool check = true;

	switch (ranCardinal)
	{
		case 0:
		{
			x1 = x[ranSite];
			y1 = (y[ranSite] + 1) % maxYLoc;
			break;
		}
		case 1:
		{
			x1 = x[ranSite];
			y1 = (y[ranSite] - 1 + maxYLoc) % maxYLoc;
			break;
		}
		case 2:
		{
			x1 = (x[ranSite] + 1) % maxXLoc;
			y1 = y[ranSite];
			break;
		}
		case 3:
		{
			x1 = (x[ranSite] - 1 + maxXLoc) % maxXLoc;
			y1 = y[ranSite];
			break;
		}
		default: std::cout << "Error checking if the current particles newly generated cell is already occupied by the current particle." << std::endl; break;
	}
	i = 0;
	while (check)
	{
		if ((x1 == x[i]) && (y1 == y[i]))
			check = false;
		i++;
		if (i == psz) break;
	}
	return check;
}

//Provides the method calls and logic to move the simulation forward a time step.
void dust_list::moveStep(int ** &updateWorld)
{

	refWorld = updateWorld;

	int myXStep, myYStep, curx, cury, stickx, sticky;
	int ptclIndx, id;
	bool xneg, check, stuck, pendingMerge, hasMoved, filter;
	int stuckCount, ptclsHandled, toRemove, hugeCount;
	int ptclSize, collGrainID;
	int numSplit;

	//Debug lines
/*	//std::cout << "Particles moving in order: " << std::endl;
	//for(unsigned int i = 0; i < order.size(); i++)
	//	std::cout << order[i] << " " << myDustList[order[i]].getID() << std::endl;
//	FILE * pFile = 0;
//	if (pFile == 0)
//		pFile = fopen("dustfile.txt", "a");   // OUTPUT: dust, size, y-position, y-step, x-localSp, y-localSp


//	FILE * pFileMoving = 0;
//	if (pFileMoving == 0)
//		pFileMoving = fopen("dustfileMoving.txt", "a");  // OUTPUT: "still moving" ptles

//	FILE * pFileStuck = 0;
//	if (pFileStuck == 0)
//		pFileStuck = fopen("dustfileStuck.txt", "a");  // OUTPUT: "still moving" ptles

//	FILE * pFilePill = 0;
//	if (pFilePill == 0)
//		pFilePill = fopen("dustfilePillCount.txt", "a");  // OUTPUT: # ptles in each  pillbox

//	FILE * pFileMerge = 0;
//	if (pFileMerge == 0)
//		pFileMerge = fopen("dustMergers.txt", "a");  // OUTPUT: # ptles in each  pillbox */

	std::vector<int> rem_dust(0);  // define a vector "remnant moving dust" in the chamber
	resetPBCounts();
	resetPotentialBlock();
	resetHandled();
	clear_width_dstr();
	clear_size_dstr();
	stuckCount = 0;
	ptclsHandled = 0;
	toRemove = 0;
	hugeCount = 0;
	numSplit = 0;

	std::vector<int> order = createRandomOrder();
	std::vector<int>::iterator iter;

	std::vector < std::vector < int > > pBoxes = getPillBoxes();
	for(int i = 0; i < myTotal; i++)
	{
		if (!myDustList[i].getFilter())
		{
			myXStep = (int)((((maxXVel + myDustList[i].getColXVel())* myGenerator->Ran()) - (maxXVel / 2.0)) / (myDustList[i].getSize()));
			myYStep = (int)(((maxYVel + myDustList[i].getColYVel())* myGenerator->Ran()) / (myDustList[i].getSize()));
			myDustList[i].setMaxXStep(myXStep);
			myDustList[i].setMaxYStep(myYStep);
			myDustList[i].setMoved(false);
			myDustList[i].setPrevPB(myDustList[i].getCurPB());
			myDustList[i].setColXVel(0);
			myDustList[i].setColYVel(0);
		}
	}

	for (iter = order.begin(); iter != order.end(); ++iter)
	{
		id = *iter;
		ptclIndx = getVecLocByID(id);
		if(ptclIndx == -1)
			std::cout << "Error: 1, starting" << std::endl;
		stuck = myDustList[ptclIndx].getStuck();
		hasMoved = myDustList[ptclIndx].checkMoved();
		pendingMerge = myDustList[ptclIndx].checkMerge();
		ptclSize = myDustList[ptclIndx].getSize();
		filter = myDustList[ptclIndx].getFilter();

		//Counts particles too large to move
		if (ptclSize > maxXVel/2 &&  ptclSize > maxYVel && !filter)
		{
			hugeCount++;
			curx = 0;
			cury = 0;
			continue;
			//std::cout << "Particle " << id << " is too large." << std::endl;
		}

		//Debug lines
/*		std::cout << "Particle " << id << " is move number " << ptclsHandled+1 << " at element " << ptclIndx << " of " << myTotal << "." << std::endl;
		std::cout << "Particle " << id << " has size " << ptclSize << "." << std::endl;
		std::cout << "Particle " << id << " has initial core location " << myDustList[ptclIndx].getXatc(0) << "," << myDustList[ptclIndx].getYatc(0) << "." << std::endl;
		std::cout << "Particle " << id << " has max X vel " << maxXVel/2 << " and max Y vel " << maxYVel << "." << std::endl;

		if(hasMoved)
			std::cout << "Has already moved." << std::endl;
		else
			std::cout << "Has not yet moved." << std::endl;

		if(pendingMerge)
			std::cout << "Is pending a merge with another particle." << std::endl;
		else
			std::cout << "Is not pending a merge with another particle." << std::endl;
		if(stuck)
			std::cout << "Is stuck to the filter" << std::endl;
		else
			std::cout << "Is not stuck to the filter." << std::endl;
*/

		//If the current particle was already flagged as stuck it will stay stuck, don't bother checking it's movement. Eventually add release probability
		else if (filter)
		{
			ptclsHandled++;
			curx = 0;
			cury = 0;
//			std::cout << id << " is a filter, and has been handled." << std::endl;
			continue;
		}
		else if (stuck)
		{
			stuckCount++;
			ptclsHandled++;
			curx = 0;
			cury = 0;
//			std::cout << "Stuck ptcl has been handled." << std::endl;
		}
		else if (pendingMerge)
		{
			ptclsHandled++;
			curx = 0;
			cury = 0;
//			std::cout << "Particle pending merge has been handled." << std::endl;
		}
		else if (!stuck && !hasMoved && !pendingMerge && !filter)
		{
			//Randomly generates the 'max' movement range for the given particle in the current time step for both x and y coordinates.
			//X being horizontal currently allows for dust to drift left or right, so the maximum horizontal movement is really half of the input 'xSpeed'
			//Effective max movement is weighted by the size of the particle, larger particles move slower.
			myXStep = myDustList[ptclIndx].getMaxXStep();
			myYStep = myDustList[ptclIndx].getMaxYStep();

			update_width_dstr(myDustList[ptclIndx].calculateWidth());
			update_size_dstr(myDustList[ptclIndx].getSize());
			//Debug line
//			std::cout << "Particle " << id << " has ran X vel " << myXStep << " and ran Y vel " << myYStep << "." << std::endl;

			//Set a boolean flag for if the x movement is left (-) or right (+)
			xneg = (myXStep < 0) ? true : false;
			curx = 0;
			cury = 0;
			stickx = 0;
			sticky = 0;
			//Checks if the location that myXStep and myYStep would move the particle into is occupied,
			//if not it finds the furthest movement possible along that vector.

			//Debug line
			//std::cout << "Checking if ptcl " << id << " can move with vel vector " << curx << "," << cury << "." << std::endl;

			check = canMakeMove(curx, cury, id);
			/*if(!check)
				std::cout << "Initial move check for " << id << " failed." << std::endl;
			else
				std::cout << "Initial move check for " << id << " succeeded." << std::endl;			*/
			while (check)
			{
				if ((curx == myXStep) && (cury == myYStep))
					break;
				if (cury < myYStep)
				{
					cury++;
					check = canMakeMove(curx, cury, id);
					if (!check)
					{
						sticky = 1;
						break;
					}
				}
				if (std::abs(curx) < std::abs(myXStep))
				{
					xneg ? curx-- : curx++;
					check = canMakeMove(curx, cury, id);
					if (!check)
					{
						xneg ? stickx = -1 : stickx = 1;
						break;
					}
				}
			}

			//Once a conflict has been found the x and y movement values are decrimented such that a viable location is held by curx and cury.
			//This might always just be -1 to each, will need to debug.
			if (!check)
			{
//				std::cout << "Particle " << id << " has come across a barrier durings it's random move vector." << std::endl;

				//If it's stopped by something in the y direction.
				if (sticky != 0)
					--cury;
				//If it's stopped by something in the x direction.
				if (stickx != 0)
					xneg ? curx++ : curx--;
//				std::cout << "Particle " << id << " has final X vel " << curx << " and final Y vel " << cury << "." << std::endl;
			}
				//After the final destination is chosen we go through all of the particles originally occupied pixels (since it's 2D) and remove them
			for (int i = 0; i < ptclSize; ++i)
				refWorld[myDustList[ptclIndx].getYatc(i)][myDustList[ptclIndx].getXatc(i)] = -1;

				//Move the core of the particle in question
			myDustList[ptclIndx].moveStep(curx, cury); // moveStep(int x, int y) is defined at dustgrain.cpp

				//Re-extend that particle in space
			for (int j = 0; j < ptclSize; ++j)
				refWorld[myDustList[ptclIndx].getYatc(j)][myDustList[ptclIndx].getXatc(j)] = id;
//				std::cout << "Particle " << id << " moved successfully to a new core location " << myDustList[ptclIndx].getXatc(0) << "," << myDustList[ptclIndx].getYatc(0) << "." << std::endl;

			if (stickx != 0 || sticky != 0)
			{
				collGrainID = getCollidingGrain(stickx, sticky, ptclIndx);
				if (collGrainID == -1)
				{
					std::cout << "An error occured: phantom grain collision." << std::endl;
				}
				else if (collGrainID == id)
				{
					std::cout << "An error occured: grain self collision." << std::endl;
				}
				else
				{
					int collGrainIndx = getVecLocByID(collGrainID);
					if (ptclIndx == -1)
						std::cout << "Error: 2, merging" << std::endl;
//					std::cout << "Particle id " << id << " and " << collGrainID << " have collided at..." << std::endl;
					if (myDustList[collGrainIndx].getFilter() && enableSticking)
					{
						myDustList[ptclIndx].setStuck(true);
						stuckCount++;
						myDustList[ptclIndx].setMoved(true);
					}
					//Comment out the else statement below in order to turn off dust-to-dust merging
					else if (enableMerging)
					{
						//We need to properly account for a particle thats too large to move merging.

						int old_size1, old_size2, new_size; // replaceX, replaceY;
						old_size1 = myDustList[ptclIndx].getSize();
						old_size2 = myDustList[collGrainIndx].getSize();
					/*	if (old_size1 > old_size2)
						{
							mergedID = id;
							myDustList[ptclIndx].setMoved(true);
							myDustList[collGrainIndx].setMerge(true);
						}
						else
						{
							mergedID = collGrainID;
							myDustList[collGrainIndx].setMoved(true);
							myDustList[ptclIndx].setMerge(true);
						}					*/

						mergeGrain_g2_to_g1(ptclIndx, collGrainIndx);
						new_size = myDustList[ptclIndx].getSize();
						update_dstr_merge(old_size1, old_size2, new_size);
						//Updates world so that old particles location get id of merged particle
						for (int c = 0; c < new_size; c++)
						{
							int replaceX = myDustList[ptclIndx].getXatc(c);
							int replaceY = myDustList[ptclIndx].getYatc(c);
							refWorld[replaceY][replaceX] = id;
						}
						toRemove++;
					}

				}
			}
			//Sets 'previous' pillbox the particle was in equal to the current pillbox it's in, as the current value has not yet been updated in this iteration.
			//TODO: Check if I can clean this all up to make it just one call for each non-stuck dust particle, something like checkContainer(xloc, yloc) where this function will save previous value, update to new value and call setPotentialBlock() as needed.
			for (unsigned int i = 0; i < pBoxes.size(); i++)
			{
				//If the particle's core is in a pillbox, add to the pillbox counter and set the dusts current pillbox then break
				//and set flag for outside any pillboxes
				if ((myDustList[ptclIndx].getXatc(0) > pBoxes[i][0]) && (myDustList[ptclIndx].getXatc(0) < pBoxes[i][1])
					&& (myDustList[ptclIndx].getYatc(0) > pBoxes[i][2]) && (myDustList[ptclIndx].getYatc(0) < pBoxes[i][3]))
				{
					//Add one to the number of particles caught by filter 'i'
					incrimentPBCounts(i);
					myDustList[ptclIndx].setCurPB(i);
					break;
				}
				else
				{
					myDustList[ptclIndx].setCurPB(-1);
					if (myDustList[ptclIndx].getPrevPB() == (int)i)
					{
						setPotentialBlock(i, false);
					}
				}
			}
			// TRACKING ONLY MOVING PTLES: => (timeCount, # of moving particles)
			// Is going to be writen in dustfileMoving.txt
			if (curx != 0 || cury != 0)
				rem_dust.push_back(id);           // Add item at end of "I am still moving" vector
			//Close parens for the stuck if statement, comment it out in the line below to stop sticking
		}//End of if(!stuck)

	//pFile=dustfile.txt
		//Commented out for optimization
	//Outputs the particles ID, size, x and y location of its core, x and y displacement during this timestep
	//fprintf(pFile, "%d %d %d %d %d %d \n",
	//id, myDustList[ptclIndx].getSize(), myDustList[ptclIndx].getXatc(0), myDustList[ptclIndx].getYatc(0), curx, cury);

	ptclsHandled++;
	}   // after done attempting move for all grains the while loop is closed

	//Counts number of ptcls to remove
	//std::cout << toRemove << " ptcles are to be removed." << std::endl;
	//Removes merged ptcls from dustList
	int oldTotal, removed;
	if (enableMerging)
	{
		oldTotal = myTotal;
		for (int r = 0; r < toRemove; r++)
		{
			removeMergedGrain();
		}
		removed = oldTotal - myTotal;
	}
	else
	{
		oldTotal = 0;
		removed = 0;
	}

	//Comment out the nest two for statements below to turn off large dust grain splitting
	if (enableSplitting)
	{
		for (int i = 0; i < myTotal; i++)
		{
			id = myDustList[i].getID();
			ptclIndx = i;
			filter = myDustList[i].getFilter();
			if (ptclIndx == -1)
				std::cout << "Error: 4, splitting" << std::endl;
			ptclSize = myDustList[i].getSize();
			if (ptclSize > maxXVel / 2 && ptclSize > maxYVel && !filter)
			{
//				std::cout << "Attempting to split up particle " << id << "." << std::endl;
				//Tell the world that the space occupied by the grain to be split is empty
				for (int c = 0; c < ptclSize; c++)
				{
					int replaceX = myDustList[ptclIndx].getXatc(c);
					int replaceY = myDustList[ptclIndx].getYatc(c);
					refWorld[replaceY][replaceX] = -1;
				}
				//Shrink existing grain and save new grain to be added later
				newlySplitGrain = attemptBreakUp(ptclIndx);
				//Tell the world which cells the original particle still occupies
				ptclSize = myDustList[ptclIndx].getSize();
				for (int c = 0; c < ptclSize; c++)
				{
					int replaceX = myDustList[ptclIndx].getXatc(c);
					int replaceY = myDustList[ptclIndx].getYatc(c);
					refWorld[replaceY][replaceX] = id;
				}
				grainsToAdd.push_back(newlySplitGrain);
				numSplit++;
			}
		}

		for (unsigned int i = 0; i < grainsToAdd.size(); i++)
		{
			myDustList.push_back(grainsToAdd[i]);
			ptclSize = myDustList[myTotal].getSize();
			id = myDustList[myTotal].getID();
			for (int c = 0; c < ptclSize; c++)
			{
				int replaceX = myDustList[myTotal].getXatc(c);
				int replaceY = myDustList[myTotal].getYatc(c);
				refWorld[replaceY][replaceX] = id;
			}
			setNewTotal();
		}
		grainsToAdd.clear();
	}
	//	std::cout << removed << " ptcles have been removed." << std::endl;

	//Sets 'previous' pillbox the particle was in equal to the current pillbox it's in, as the current value has not yet been updated in this iteration.
	//	scanPBoxes();
	std::vector < bool > poreFilled = checkPoreFilled();

	for (unsigned int i = 0; i < pillBoxes.size(); i++)
	{
		if (!poreFilled[i])
		{
			//Tells the pore that the particle just exited that it can not be clogged as of the current timestep.
			setPotentialBlock(i, false);
		}
	}

	//Combines the information about the pore's line density and if it has had a zero flux for some set time limit, allowing us to roughly determine whether a pore is actually blocked.
	checkBlocked();

	//Writes the number of particles inside each pillbox to a text file.
/*	for (unsigned int i = 0; i < pillBoxes.size(); i++)
	{
		fprintf(pFilePill, "%d ", // pFile=dustfilePillCount.txt
				pBCounts[i]);
	}
	fprintf(pFilePill, "\n");
					*/
	//Writes out the number of particles stuck to the filters to a text file.
	//fprintf(pFileStuck, "%i \n", stuckCount); //removed for optimizing

	unsigned int rem_dustS = rem_dust.size();  // the size of this vector is the number of "still moving" ptles
	//Writes out number of moving, stuck, too large to move, merged, and total


	//Edited so that dust sizes are only printed out once, for cases where merging and splitting are disabled
	int curTime = getTimeSteps();
	if (curTime % 100 == 0)
	{
//		std::vector <double> dustStats(0);
		std::vector <double> sizeStats(0);
		std::vector <double> widthStats(0);

		sizeStats = calc_stats(sizeDist);
		widthStats = calc_stats(dustWidth);
//		dustStats = calc_stats(dustDist);

		FILE * pFileD = 0;
		std::string dustDistStatFile = procOutputFolder + "/dustDistStats.txt";
		if (pFileD == 0)
			pFileD = fopen(dustDistStatFile.c_str(), "a");  // OUTPUT: # ptles sizes
		fprintf(pFileD, "%i %f %f %f %f \n", curTime, widthStats[0], widthStats[1], sizeStats[0], sizeStats[1]);
		fclose(pFileD);

		FILE * pFileW = 0;
		std::string dustWidthListFile = procOutputFolder + "/dustWidthList.txt";
		if (pFileW == 0)
			pFileW = fopen(dustWidthListFile.c_str(), "a");
		fprintf(pFileW, "%i", curTime);
		for(unsigned int i = 0; i < dustWidth.size(); i++)
		{
																								int width = dustWidth[i][0];
																								int wCount = dustWidth[i][1];
			fprintf(pFileW, " %i %i", width, wCount);

		}
		fprintf(pFileW, "\n");
		fclose(pFileW);

		FILE * pFileS = 0;
		std::string dustSizeListFile = procOutputFolder + "/dustSizeList.txt";
		if (pFileS == 0)
			pFileS = fopen(dustSizeListFile.c_str(), "a");
		fprintf(pFileS, "%i", curTime);
		for(unsigned int i = 0; i < sizeDist.size(); i++)
		{
																								int size = sizeDist[i][0];
																								int sCount = sizeDist[i][1];
			fprintf(pFileS, " %i %i ", size, sCount);
		}
		fprintf(pFileS, "\n");
		fclose(pFileS);

		FILE * cFile = 0;
		std::string dustCountFile = procOutputFolder + "/dustCount.txt";
		if (cFile == 0)
			cFile = fopen(dustCountFile.c_str(), "a");   // OUTPUT: dust, size, y-position, y-step, x-localSp, y-localSp
		fprintf(cFile, "%i %i %i %i %i %i \n",(int) rem_dustS, stuckCount, hugeCount, numSplit, removed, myTotal);
		fclose(cFile);
	}
	//fprintf(cFile, "%d \n", myTotal); //Totalgrains remaining removed for optimizing
	//fprintf(pFileMerge, "%d \n", removed);

	//Move's the dust_list's time counter one timestep forward.
	incrimentTimeStep();

	rem_dust.clear();

//	fclose(pFilePill);


//	fclose(pFileMerge);
//	fclose(pFile);
//	fclose(pFileMoving);
//	fclose(pFileStuck);

	//////////////////////////////////////////////////////////
}
void dust_list::calculatePostCollisionVelocities(dust_grain movingPtcl, dust_grain staticPtcl, int movActXVel, int movActYVel)
{
//	int statXVelBefore = staticPtcl.getPrevXVel();
//	int statYVelBefore = staticPtcl.getPrevYVel();
//	int statXMomBefore = staticPtcl.getSize()*statXVelBefore;
//	int statYMomBefore = staticPtcl.getSize()*statYVelBefore;
	int movXVelBefore = movingPtcl.getMaxXStep();
	int movYVelBefore = movingPtcl.getMaxYStep();
//	int movXMomBefore = movingPtcl.getSize()*movXVelBefore;
//	int movYMomBefore = movingPtcl.getSize()*movYVelBefore;
//	int netXMom = movXMomBefore + statXMomBefore;
//	int netYMom = movYMomBefore + statYMomBefore;
//	int statXMomAfter, statYMomAfter, movXMomAfter, movYMomAfter, statXVelAfter, statYVelAfter, movXVelAfter, movYVelAfter;

	int xImpulse = (movXVelBefore - movActXVel)*movingPtcl.getSize();
	int yImpulse = (movYVelBefore - movActYVel)*movingPtcl.getSize();

	staticPtcl.setColXVel(xImpulse/staticPtcl.getSize());
	staticPtcl.setColYVel(yImpulse/staticPtcl.getSize());
	
	//Check if it hit statPtcl from behind or the side


}

void dust_list::setFunctionality(bool splitting, bool sticking, bool merging)
{
	enableSplitting = splitting;
	enableSticking = sticking;
	enableMerging = merging;
}

//Can save time and make things more efficient by making the dust's location vectors a tuple, so they can be manipulated simultaneously and reordered.
dust_grain dust_list::attemptBreakUp(int grain)
{
	dust_grain dCopy = myDustList[grain];
	int mySize = dCopy.getSize();
	std::vector < int > firstShardX;
	std::vector < int > firstShardY;
	std::vector < int > secondShardX;
	std::vector < int > secondShardY;
	int hX = 0;
	int hY = 0;
	int lX = maxXLoc * 2;
	int lY = maxYLoc * 2;
	int xPos;
	int yPos;
	bool boundary;
	std::vector < int > copyMyX = std::vector < int >(mySize);
	std::vector < int > copyMyY = std::vector < int >(mySize);
	//Create a copy of continuous dust occupied cells to account for a grain being partially across a boundary
	for (int i = 0; i < mySize; i++)
	{
		xPos = dCopy.getXatc(i);
		copyMyX[i] = xPos;
		if (xPos < lX)
			lX = xPos;
		if (xPos > hX)
			hX = xPos;

		yPos = dCopy.getYatc(i);
		copyMyY[i] = yPos;
		if (yPos < lY)
			lY = yPos;
		if (yPos > hY)
			hY = yPos;

		//Debug statement to make sure that particles at the boundaries are handled properly
		if (xPos <= 1 || yPos <= 1)
			boundary = true;
		else
			boundary = false;
	}
	if (boundary)
		std::cout << "Particle " << dCopy.getID() << " is lying on a boundary during split." << std::endl;

	bool xEdge = false;
	int xDiff = hX - lX;
	if (xDiff >= mySize)
	{
		xEdge = true;
		for (int i = 0; i < mySize; i++)
		{
			if (copyMyX[i] < (mySize-1))
			{
				copyMyX[i] = copyMyX[i] + maxXLoc;
			}
		}
	}
	bool yEdge = false;
	int yDiff = hY - lY;
	if (yDiff >= mySize)
	{
		yEdge = true;
		for (int i = 0; i < mySize; i++)
		{
			if (copyMyY[i] < (mySize-1))
			{
				copyMyY[i] = copyMyY[i] + maxYLoc;
			}
		}
	}
	if (xEdge)
	{
		lX = maxXLoc * 2;
		hX = 0;
		for (int i = 0; i < mySize; i++)
		{
			xPos = copyMyX[i];
			if (xPos < lX)
				lX = xPos;
			if (xPos > hX)
				hX = xPos;
		}
	}
	if(yEdge)
	{
		lY = maxYLoc * 2;
		hY = 0;
		for (int i = 0; i < mySize; i++)
		{
			yPos = copyMyY[i];
			if (yPos < lY)
				lY = yPos;
			if (yPos > hY)
				hY = yPos;
		}
	}
	int temp;
	int countA = 0;
	int countB = 0;

	if (xDiff > yDiff)
	{
		for (unsigned int i = 0; i < copyMyX.size(); i++)
		{
			if (copyMyX[i] == lX && countA == 0)
			{
				firstShardX.push_back(copyMyX[i]);
				firstShardY.push_back(copyMyY[i]);
				temp = i;
				countA++;
				break;
			}
		}
		copyMyX.erase(copyMyX.begin() + temp);
		copyMyY.erase(copyMyY.begin() + temp);

		for (unsigned int i = 0; i < copyMyX.size(); i++)
		{
			if (copyMyX[i] == hX && countB == 0)
			{
				secondShardX.push_back(copyMyX[i]);
				secondShardY.push_back(copyMyY[i]);
				temp = i;
				countB++;
				break;
			}
		}
		copyMyX.erase(copyMyX.begin() + temp);
		copyMyY.erase(copyMyY.begin() + temp);
	}
	else
	{
		for (unsigned int i = 0; i < copyMyX.size(); i++)
		{
			if (copyMyY[i] == lY && countA == 0)
			{
				firstShardX.push_back(copyMyX[i]);
				firstShardY.push_back(copyMyY[i]);
				temp = i;
				countA++;
				break;
			}
		}
		copyMyX.erase(copyMyX.begin() + temp);
		copyMyY.erase(copyMyY.begin() + temp);

		for (unsigned int i = 0; i < copyMyX.size(); i++)
		{
			if (copyMyY[i] == hY && countB == 0)
			{
				secondShardX.push_back(copyMyX[i]);
				secondShardY.push_back(copyMyY[i]);
				temp = i;
				countB++;
				break;
			}
		}
		copyMyX.erase(copyMyX.begin() + temp);
		copyMyY.erase(copyMyY.begin() + temp);
	}

	int tries;
	while ((countA + countB) < mySize)
	{
		if ((countA <= countB))
		{
			core:
			tries = 0;
			for (unsigned int i = 0; i < copyMyX.size(); i++)
			{
				for (unsigned int j = 0; j < firstShardX.size(); j++)
				{
					if (((copyMyX[i] == firstShardX[j]) && ((copyMyY[i] == firstShardY[j] + 1) || (copyMyY[i] == firstShardY[j] - 1)))
						|| ((copyMyY[i] == firstShardY[j]) && ((copyMyX[i] == firstShardX[j] + 1) || (copyMyX[i] == firstShardX[j] - 1))))
					{
						firstShardX.push_back(copyMyX[i]);
						firstShardY.push_back(copyMyY[i]);
						temp = i;
						countA++;
						goto next;
					}
				}
				tries++;
				if (tries == (int)copyMyX.size())
				{
					goto shard;
				}
			}
		}
		else if ((countA > countB))
		{
			shard:
			tries = 0;
			for (unsigned int i = 0; i < copyMyX.size(); i++)
			{
				for (unsigned int j = 0; j < secondShardX.size(); j++)
				{
					if (((copyMyX[i] == secondShardX[j]) && ((copyMyY[i] == secondShardY[j] + 1) || (copyMyY[i] == secondShardY[j] - 1)))
						|| ((copyMyY[i] == secondShardY[j]) && ((copyMyX[i] == secondShardX[j] + 1) || (copyMyX[i] == secondShardX[j] - 1))))
					{
						secondShardX.push_back(copyMyX[i]);
						secondShardY.push_back(copyMyY[i]);
						temp = i;
						countB++;
						goto next;
					}
				}
				tries++;
				if (tries == (int)copyMyX.size())
				{
					goto core;
				}
			}

		}
	next:
		copyMyX.erase(copyMyX.begin() + temp);
		copyMyY.erase(copyMyY.begin() + temp);
	}

	std::vector <dust_grain> splitDust;

	for (unsigned int i = 0; i < firstShardX.size(); i++)
	{
		firstShardX[i] = firstShardX[i] % maxXLoc;
		firstShardY[i] = firstShardY[i] % maxYLoc;
	}
	for (unsigned int i = 0; i < secondShardX.size(); i++)
	{
		secondShardX[i] = secondShardX[i] % maxXLoc;
		secondShardY[i] = secondShardY[i] % maxYLoc;
	}

	dCopy.growGrain(firstShardX, firstShardY);
	dCopy.setStuck(false);
	dCopy.setMoved(true);
	myDustList[grain] = dCopy;

	dust_grain shard = dust_grain(secondShardX, secondShardY, secondShardX.size(), newUniqueID());

	shard.setStuck(false);
	shard.setMoved(true);
	shard.setPrevPB(dCopy.getPrevPB());
	shard.setCurPB(dCopy.getCurPB());
	shard.setFilter(false);

	update_dstr_split(mySize, dCopy.getSize(), shard.getSize());

	return shard;

}
void dust_list::incrimentTimeStep()
{
	numTimeSteps++;
}

/*  getCollidingGrain was from original version handed to me at project start. Not yet implimented, but has been tweaked to reflect other major changes to the code base.
	Method to merge two regular dust grains/particles
	Not yet tested or implimented
	Need to allow for collisions in the canMakeMove method	*/
int dust_list::getCollidingGrain(int xmove, int ymove, int grain_self) //KM
{
	int collgrain, id, j, x, y;
	dust_grain temp;
	temp = getGrainByVecLoc(grain_self);
	j = temp.getSize();
	for (int c = 0; c < j; c++)
	{
		x = (temp.getXatc(c) + xmove + maxXLoc) % maxXLoc;
		y = (temp.getYatc(c) + ymove + maxYLoc) % maxYLoc;

		id = refWorld[y][x];
		if ((id != -1) && (id != temp.getID()))
		{
			collgrain = id;
			return collgrain;
		}
	}
	return -1;
}

void dust_list::setNewTotal()
{
	myTotal = myDustList.size();
}

std::vector < bool > dust_list::checkPoreFilled()
{
	int xloc, yloc, tally;
	std::vector < bool > outputBool;
	std::vector< std::vector < int > > pBoxes = getPillBoxes();
	int loops = pBoxes.size();
	for (int i = 0; i < loops; i++)
	{
		tally = 0;
		xloc = pBoxes[i][0];
		yloc = pBoxes[i][2];
		int diff = pBoxes[i][1] - pBoxes[i][0];
		for (int j = 0; j < diff; j++)
		{
			if (isOpen(xloc, yloc))
			{
				tally++;
			}
			xloc = (int)(xloc + 1) % maxXLoc;
		}
		int maxEmpty = int(diff*0.7);
		if (tally < maxEmpty)
			outputBool.push_back(true);
		else
			outputBool.push_back(false);
	}
	return outputBool;
}

//Checks if a dust particle can move to the space indicated by a given x and y velocity.
bool dust_list::canMakeMove(int xmove, int ymove, int grainID)
{
	bool canMakeIt = true;
	int indx = getVecLocByID(grainID);
	if (indx == -1)
		std::cout << "Error: 3, moving" << std::endl;
	int size = myDustList[indx].getSize();
//	std::cout << "Checking if dust id " << grainID << " can move with vector: " << xmove << "," << ymove << "." << std::endl;
	for (int c = 0; c < size; c++)
	{
		/* The following controls overlaping.
		If !false, particles are transparent each other. Also the filter
			*/

		//The modulus calculations are for the wrap around boundary condition, but the addition should only affect the x calculation since it could have a negative value without said addition.
		if (!isOpen(((myDustList[indx].getXatc(c) + xmove + maxXLoc) % maxXLoc), ((myDustList[indx].getYatc(c) + ymove + maxYLoc) % maxYLoc), grainID))
			canMakeIt = false;

	}

	return canMakeIt;

	/*Alternative: Check xmove and ymove sepaarately
	if(!isOpen(( (myDustList[grainNumber].getXatc(c)+maxXLoc)%maxXLoc), ((myDustList[grainNumber].getYatc(c)+ymove+maxYLoc)%maxYLoc),grainNumber))
							return false;
	if(!isOpen(( (myDustList[grainNumber].getXatc(c)+xmove+maxXLoc)%maxXLoc), ((myDustList[grainNumber].getYatc(c)+maxYLoc)%maxYLoc),grainNumber))
							return false;	*/

}

//Actually merges dust particles to filter, removing them from the dust_list and adding a new larger particle at the location of the merger.
//Currently not used, there is instead a flag that gets set in moveStep if the center of a particle impacts with the filter.
//TODO: Impliment dust sticking via this method. It would likely be more efficient but would make it very hard to impliment un-sticking.
/*dust_grain dust_list::mergeGrain_to_filter(int g1, int g2) //KM
{
	dust_grain temp, g1grain, g2grain;
	int g1sze, g2sze, totsze, c;

	g1grain = myDustList[g1];
	g2grain = myDustList[g2];
	g1sze = g1grain.getSize();
	g2sze = g2grain.getSize();
	std::cout << "Sizes are " << g1sze << "\t" << g2sze << std::endl;
	totsze = g1sze + g2sze;
	std::vector<int> x(totsze);
	std::vector<int> y(totsze);
	for (c = 0; c < g1sze; c++)
	{
		x[c] = g1grain.getXatc(c);
		y[c] = g1grain.getYatc(c);
		std::cout << x[c] << "\t" << y[c] << std::endl;
	}
	std::cout << std::endl;
	for (c = 0; c < g2sze; c++)
	{
		x[c + g1sze] = g2grain.getXatc(c);
		y[c + g1sze] = g2grain.getYatc(c);
		std::cout << x[c] << "\t" << y[c] << std::endl;
	}
	temp = dust_grain(x, y, totsze);
	return temp;
}*/
void dust_list::setProcOutputFolder(std::string dirName)
{
	procOutputFolder = dirName;
}

//Actually merges two dust particles, removing them from the dust_list and adding a new larger particle at the location of the merger.
//This must be called after the dust to dust impact calculation is performed.
void dust_list::mergeGrain_g2_to_g1(int g1Indx, int g2Indx) //KM
{
	dust_grain g1grain, g2grain;
	int g1sze, g2sze, totsze, c;
	g1grain = myDustList[g1Indx];
	g2grain = myDustList[g2Indx];
	g1sze = g1grain.getSize();
	g2sze = g2grain.getSize();
//	std::cout << "Sizes of particles to merge are " << g1sze << " and " << g2sze << "." << std::endl;
	totsze = g1sze + g2sze;
	std::vector<int> x(totsze);
	std::vector<int> y(totsze);
	for (c = 0; c < g1sze; c++)
	{
		x[c] = g1grain.getXatc(c);
		y[c] = g1grain.getYatc(c);
//		std::cout << x[c] << "\t" << y[c] << std::endl;
	}
//	std::cout << std::endl;
	for (c = 0; c < g2sze; c++)
	{
		x[c + g1sze] = g2grain.getXatc(c);
		y[c + g1sze] = g2grain.getYatc(c);
//		std::cout << x[c + g1sze] << "\t" << y[c + g1sze] << std::endl;
	}
	myDustList[g1Indx].growGrain(x, y);
	myDustList[g1Indx].setSize();
	myDustList[g1Indx].setMoved(true);
	myDustList[g2Indx].setMerge(true);
	myDustList[g2Indx].clearGrain();

}

/*	//Method which allows for a harmonic shaking of the filter, will need to step through and clean up this routine.
//One more function that isn't yet implimented in this version... will also need heavy editing to get working.
void dust_list::addGrainSk(int filterGap, int filterWidth, int filterLength, float timeF, int ** &updateWorld)
{

	//refWorld = updateWorld;  //talves no necesito esto "perhaps they do not need this"

	std::cout << "Adding Filter" << std::endl;
	int Chunk = filterGap + filterWidth;
	int numChunk = maxXLoc / Chunk;
	int RemChunk = maxXLoc % Chunk;
	int cells = numChunk * filterWidth * filterLength;
	int id =
	if (RemChunk > 0)
	{
		cells += std::max(filterWidth, RemChunk - filterGap)*filterLength;

	}

	std::vector<int> x = std::vector<int>(cells);
	std::vector<int> y = std::vector<int>(cells);

	//"timeF" must be a fuction of time.
	//2.1, to down, 1.9 to top

	int start = int(maxYLoc / (2.5 + (0.18*timeF)));
	int fG, fW, count;
	bool on = true;

	fG = fW = 0;
	count = 0;

	for (int i = 0; i < filterLength; ++i)
	{
		on = true;
		fG = fW = 0;
		for (int j = 0; j < maxXLoc; ++j)
		{
			if (on)
			{
				fW++;
				y[count] = i + start;
				x[count] = j;
				++count;

				if (fW == filterWidth)
				{
					on = false;
					fW = 0;
				}
			}
			else
			{
				fG++;

				if (fG == filterGap)
				{
					on = true;
					fG = 0;
				}
			}
		}
	}


	//IncreaseListbyOne();  // Here we don't actually change the size of the list, we're just shifting the filter back and forth.

	dust_grain temp = dust_grain(x, y, cells, id);
	myDustList[myTotal - 1] = temp;

	for (int c = 0; c < cells; ++c)
		refWorld[y[c]][x[c]] = id;

}	*/

//Creates the first filter based on input parameters and inserts it in the simulation.
void dust_list::addGrain(int filterGap, int filterWidth, int filterLength)
{
	std::cout << "Adding Filter..." << std::endl;
	int Chunk = filterGap + filterWidth;
	int numChunk = maxXLoc / Chunk;
	int RemChunk = maxXLoc % Chunk;
	int cells = numChunk * filterWidth * filterLength;
	int id = newUniqueID();
	if (RemChunk > 0)
	{
		cells += std::max(filterWidth, RemChunk - filterGap)*filterLength;
	}

	std::vector<int> x = std::vector<int>(cells);
	std::vector<int> y = std::vector<int>(cells);

	int start = maxYLoc / 2;
	int fG, fW, count;
	bool on = true;

	fG = fW = 0;
	count = 0;

	for (int i = 0; i < filterLength; ++i)
	{
		on = true;
		fG = fW = 0;
		for (int j = 0; j < maxXLoc; ++j)
		{
			if (on)
			{
				fW++;
				y[count] = i + start;
				x[count] = j;
				++count;

				if (fW == filterWidth)
				{
					on = false;
					fW = 0;
				}
			}
			else
			{
				fG++;

				if (fG == filterGap)
				{
					on = true;
					fG = 0;
				}
			}
		}
	}


	IncreaseListbyOne();
	dust_grain temp = dust_grain(x, y, cells, id);
	temp.setFilter(true);
	temp.setMaxXStep(0);
	temp.setMaxYStep(0);
	temp.setMaxXLoc(maxXLoc);
	temp.setMaxYLoc(maxYLoc);
	myDustList[myTotal - 1] = temp;

	for (int c = 0; c < cells; ++c)
		refWorld[y[c]][x[c]] = id;
	//Defines edges of the pillboxes around each pore.
	calcPillBoxes(filterGap, filterWidth, filterLength, start);
}

//Creates the second filter based on a second set of input parameters and inserts it in the simulation.
//TODO: Need to just change the code above to have it placed at a new location by the same function, maybe have the y loc determined by current size of refWorld.
void dust_list::addGrain2(int filter2Gap, int filter2Width, int filter2Length)
{
	std::cout << "Adding Filter" << std::endl;
	int Chunk = filter2Gap + filter2Width;
	int numChunk = maxXLoc / Chunk;
	int RemChunk = maxXLoc % Chunk;
	int cells = numChunk * filter2Width * filter2Length;
	int id = newUniqueID();
	if (RemChunk > 0)
	{
		cells += std::max(filter2Width, RemChunk - filter2Gap)*filter2Length;
	}


	std::vector<int> x = std::vector<int>(cells);
	std::vector<int> y = std::vector<int>(cells);


	int start = (int) 3*maxYLoc / 4;       //Vertical position of filter
	int fG, fW, count;
	bool on = true;

	fG = fW = 0;
	count = 0;

	for (int i = 0; i < filter2Length; ++i)
	{
		on = true;
		fG = fW = 0;
		for (int j = 0; j < maxXLoc; ++j)
		{

			if (on)
			{
				fW++;
				y[count] = i + start;
				x[count] = j;
				++count;

				if (fW == filter2Width)
				{
					on = false;
					fW = 0;
				}
			}
			else
			{
				fG++;

				if (fG == filter2Gap)
				{
					on = true;
					fG = 0;
				}
			}
		}
	}

	IncreaseListbyOne();
	dust_grain temp = dust_grain(x, y, cells, id);
	temp.setFilter(true);
	temp.setMaxXStep(0);
	temp.setMaxYStep(0);
	temp.setMaxXLoc(maxXLoc);
	temp.setMaxYLoc(maxYLoc);
	myDustList[myTotal - 1] = temp;

	for (int c = 0; c < cells; ++c)
		refWorld[y[c]][x[c]] = id;
	calcPillBoxes(filter2Gap, filter2Width, filter2Length, start);
}

/*	It should be a simple switch...
you can change what generator it uses... just switch the
myGenerator->Ran()
with
myGenerator->Gauss()
Only required in the first line of dust_list::addGrain	*/

double dust_list::Abs(double Nbr)
{
	//	return (Nbr >= 0) ? Nbr : -Nbr;
	if (Nbr >= 0)
		return Nbr;
	else
		return -Nbr;
}

//Randomly generataes a grain between the min and max grain size input parameters then randomly places it at an empty location in the simulation.
void dust_list::addGrain(int low, int high)
{
	//Comment out to switch between gaussian and normal particle size distributions
	//int Tsize = ((int)(1000000000*myGenerator->Ran()))%(high-low+1)+low;
	int Tsize = (int)((myGenerator->Gauss2())*(high - low)) + low;

	int ranSite, ranCardinal;
	double tRand;
	int count = 0;
	int tries = 0;
	//Sets a maximum number of times to try and extend the current particle before picking a new 'core' location.
	int maxTries = 1000;
	int id = newUniqueID();
	bool foundSpot = false;
	bool validCore = false;

	//Creates a pair of Tsize x and y values which will be the spatial locations of the current particles cells
	std::vector<int> x = std::vector<int>(Tsize);
	std::vector<int> y = std::vector<int>(Tsize);
	//We assume that the first core is invalid to get things started, if the particle finds a location to spawn we set validCore = true and break from the while to actually add the particle.
	while (!validCore)
	{
		//First generates new random initial point for the 'core' of the dust particle until an unoccupied location is found
		do
		{
			x[0] = ((int)(1000000000 * myGenerator->Ran())) % maxXLoc;
			y[0] = ((int)(1000000000 * myGenerator->Ran())) % maxYLoc;
		} while (!isOpen(x[0], y[0]));

		count = 0;
		//Now extends the dust particle in space to take up 'Tsize' cells using a normal distribution to extend one cell at a time into the cardinal neighborhood
		for (int c = 1; c < Tsize; c++)
		{
			//Continues until an empty neighboring cell is found to extend the particle into
			tries = 0;
			foundSpot = false;
			do {
				//Ran() returns a value between 0 and 1, so that the inverse will be greater than or equal to 1, but we need to make sure that we don't divide by 0
				tRand = 0;
				while (tRand == 0)
					tRand = myGenerator->Ran();

				//The 'ranSite' represents one of the cells occupied by the current particle
				ranSite = ((int)(1 / tRand)) % c;

				if (!foundSpot)
				{
					//We make sure our next random number isn't zero to avoid producing NaN by division
					tRand = 0;
					while (tRand == 0)
						tRand = myGenerator->Ran();

					//Chooses a random cardinal direction of value 0 through 3
					ranCardinal = ((int)(1 / tRand)) % 4;
					switch (ranCardinal)
					{
						//Tries to extend to one particle above the current cell 'ranSite' of the current dust particle
					case 0: {
						/*The purpose of taking the modulus of the location with the max_Loc's  (% max_Loc) is to create the wrap around boundary condition.
						A particle with y[0] == myYMax that tries to spawn a particle above it will actually then occupy y[1] >= 1)
						First logical statement simply checks if the generated cell is unoccupied by existing particles
						Second logical statement checks if the generated cell is unoccupied by the current particle, since it has not yet been added to the world it won't appear in the first check.	*/
						if ((isOpen(x[ranSite], (y[ranSite] + 1) % maxYLoc)) && (isOpenSelf(x, y, ranSite, ranCardinal, c)))
						{
							foundSpot = true;
							x[c] = x[ranSite];
							y[c] = (y[ranSite] + 1) % maxYLoc;
						}
						break;
					}
					//Tries to extend to one particle below the current cell 'ranSite' of the current dust particle
					case 1:
					{
						if ((isOpen(x[ranSite], (y[ranSite] - 1 + maxYLoc) % maxYLoc)) && (isOpenSelf(x, y, ranSite, ranCardinal, c)))
						{
							foundSpot = true;
							x[c] = x[ranSite];
							y[c] = (y[ranSite] - 1 + maxYLoc) % maxYLoc;
						}
						break;
					}
					//Tries to extend to one particle to the right of the current cell 'ranSite' of the current dust particle
					case 2: {
						if ((isOpen((x[ranSite] + 1) % maxXLoc, y[ranSite])) && (isOpenSelf(x, y, ranSite, ranCardinal, c)))
						{
							foundSpot = true;
							x[c] = (x[ranSite] + 1) % maxXLoc;
							y[c] = y[ranSite];
						}
						break;
					}
							//Tries to extend to one particle to the left of the current cell 'ranSite' of the current dust particle
					case 3: {
						if ((isOpen((x[ranSite] - 1 + maxXLoc) % maxXLoc, y[ranSite])) && (isOpenSelf(x, y, ranSite, ranCardinal, c)))
						{
							foundSpot = true;
							x[c] = (x[ranSite] - 1 + maxXLoc) % maxXLoc;
							y[c] = y[ranSite];
						}
						break;
					}
					default: std::cout << "Error extending particle " << myTotal << " in space." << std::endl; break;
					}//End of switch statement
				}//End of if(!foundSpot)
				++tries;
			//Attempts to find a new adjacent cell to place the 'c'th cell for 'maxTries' times.
			} while (!foundSpot && tries < maxTries);

			//If we stopped placing cells because we tried the max number of times we break out of the for statement and flag the 'core' as invalid and choose a new 'core' location.
			if (tries == maxTries)
			{
				validCore = false;
				break;
			}
			else
				count++;
		}//End of for (int c = 1; c < Tsize; c++)

		//If we've placed the desired number of cells for the current particle we change the flag to true and break back to and thus out of the while(!validCore)
		if (count == Tsize-1)
		{
			validCore = true;
			break;
		}
	}
	//If we've found a valid configuration we place it in the world.
	IncreaseListbyOne();
	dust_grain temp = dust_grain(x, y, Tsize, id);
	temp.setMaxXLoc(maxXLoc);
	temp.setMaxYLoc(maxYLoc);
	temp.setMaxXStep(0);
	temp.setMaxYStep(0);
	temp.setPrevXVel(0);
	temp.setPrevYVel(0);

	myDustList[myTotal-1] = temp;
	for (int c = 0; c < Tsize; ++c)
		refWorld[y[c]][x[c]] = id;
	//std::cout << "Grain id " << id << " added to world." << std::endl;

}

//Removes a specificly indexed grain from the list.
void dust_list::removeMergedGrain()
{
	int loc = -1;
	for (int c = 0; c < myTotal; c++)
	{
		if (myDustList[c].checkMerge())
		{
			loc = c;
			break;
		}
	}
	for (int indx = loc; indx < (myTotal - 1); ++indx)
		myDustList[indx] = myDustList[indx + 1];
	ShrinkListbyOne();
}

//Lowers the size of the list by one both in memory and in the int value of the total num ptcls.
void dust_list::ShrinkListbyOne()
{
	if (myTotal > 0)
	{
		myTotal--;
		myDustList.resize(myTotal);
	}

}

//Increase the size of the list by one and update myTotal.
void dust_list::IncreaseListbyOne()
{
	myTotal++;
	myDustList.resize(myTotal);
}

//Randomly selects the order in which dust_grains will perform their random movements.
std::vector<int> dust_list::createRandomOrder()
{
	std::vector<int> temp(myTotal);
	for (int i = 0; i < myTotal; ++i)
		temp[i] = getIDByVecLoc(i);

	int start = 0;
	int randomSpot;
	int tempLoc;

	while (start < myTotal)
	{
		randomSpot = (int)(1000000000 * myGenerator->Ran()) % (myTotal - start) + start;

		tempLoc = temp[start];
		temp[start] = temp[randomSpot];
		temp[randomSpot] = tempLoc;

		++start;
	}
	return temp;
}

//Should be called from addGrain functions that create filters and will create a vector of (miny, maxy, minx, maxx) co-ordinates defining the bounds of each pillbox.
void dust_list::calcPillBoxes(int filterGap, int filterWidth, int filterLength, int start)
{
	int loc = 0;
	int count = 0;

	//points[0] is x position of left corners, points[1] is x position of right corners, points[2] is y position of bottom corners, points[3] is y position of top corners
	std::vector<int> points = std::vector<int>(4);

	points[0] = points[1] = 0;
	points[2] = start - 3;
	points[3] = start + filterLength;
	//The box defined is 3 cells left, right and below the cell, but flush with the top.
	while (loc < maxXLoc)
	{
		loc = loc + filterWidth + filterGap;
		if (loc > maxXLoc)
			break;

		points[0] = points[1] + filterWidth - 3;
		points[1] = points[0] + filterGap + 6;

		dust_list::pushBackpBVectors(points);
		points[0] = points[0] + 3;
		points[1] = points[1] - 3;
		//std::cout << pillBoxes[count][0] << " " << pillBoxes[count][1] << " " << pillBoxes[count][2] << " " << pillBoxes[count][3] << " \n";

		count++;
	}
	//std::cout << "Debuggin' pillboxes." << std::endl;
}

//Lengthens vectors relevant to pore tracking/clogging, called once for each pore during their calculation/creation.
void dust_list::pushBackpBVectors(std::vector < int > corners)
{
	dust_list::pillBoxes.push_back(corners);
	dust_list::pBCounts.push_back(0);
	dust_list::poreJamTimer.push_back(0);
	dust_list::potentialBlock.push_back(true);
	dust_list::poreBlocked.push_back(false);
}

//Checks all of the pores to see if they have a certain line density of occupied cells across their width, and if they have had no flux for a set amount of time. If requirements are met the pore is flagged as clogged.
void dust_list::checkBlocked()
{
	FILE * pBlocked = 0;
	std::string poreBlockFile = procOutputFolder + "/poresBlocked.txt";
	if (pBlocked == 0)
	{
		pBlocked = fopen(poreBlockFile.c_str(), "a");   // OUTPUT: pore, timeblocked
	}

	std::vector < bool > potBlock = getPotentialBlock();
	std::vector < int > poreJamT = getPoreJamTimer();
	int time = getTimeSteps();
	int timeBlocked = 0;
	std::vector < bool > poreFilled = checkPoreFilled();
	for (unsigned int i = 0; i < potBlock.size(); i++)
	{
		if (getPoreBlocked(i) == true)
			continue;

		if (potBlock[i])
			incrimentPoreJamCounter(i);
		else
		{
			setPoreJamTimer(i, 0);
			continue;
		}

		if (poreJamT[i] > 300)
		{
			setPoreBlocked(i, true);
			timeBlocked = time - 300;
			fprintf(pBlocked, "%i %i \n", i, timeBlocked);
		}
	}
	fclose(pBlocked);
}

std::vector< std::vector< int > > dust_list::getPillBoxes()
{
	return pillBoxes;
}

int dust_list::getpBCounts(int loc)
{
	return pBCounts[loc];
}

std::vector < int > dust_list::getPoreJamTimer()
{
	return poreJamTimer;
}

bool dust_list::getPoreBlocked(int loc)
{
	return poreBlocked[loc];
}

std::vector < bool > dust_list::getPotentialBlock()
{
	return potentialBlock;
}

void dust_list::setPillBoxes(int loc, std::vector< int > pBoxes)
{
	pillBoxes[loc] = pBoxes;
}

void dust_list::setpBCounts(int loc, int pCounts)
{
	pBCounts[loc] = pCounts;
}

void dust_list::setPoreJamTimer(int loc, int jamTimer)
{
	poreJamTimer[loc] = jamTimer;
}

void dust_list::setPoreBlocked(int loc, bool verBlock)
{
	poreBlocked[loc] = verBlock;
}

void dust_list::setPotentialBlock(int loc, bool potBlock)
{
	potentialBlock[loc] = potBlock;
}

void dust_list::incrimentPBCounts(int loc)
{
	++pBCounts[loc];
}

void dust_list::incrimentPoreJamCounter(int loc)
{
	++poreJamTimer[loc];
}

void dust_list::resetPBCounts()
{
	for (unsigned int i = 0; i < pBCounts.size(); i++)
	{
		setpBCounts(i, 0);
	}
}

void dust_list::setMaxXVel(int xvel)
{
	maxXVel = xvel;
}

void dust_list::setMaxYVel(int yvel)
{
	maxYVel = yvel;	
}

int dust_list::getMaxXVel()
{
	return maxXVel;
}

int dust_list::getMaxYVel()
{
	return maxYVel;
}

void dust_list::resetPotentialBlock()
{
	for (unsigned int i = 0; i < potentialBlock.size(); i++)
	{
		setPotentialBlock(i, true);
	}
}

int dust_list::getTimeSteps()
{
	return numTimeSteps;
}

int dust_list::getMaxXLoc()
{
	return maxXLoc;
}

int dust_list::getMaxYLoc()
{
	return maxYLoc;
}

void dust_list::setMaxXLoc(int maxX)
{
	maxXLoc = maxX;
}

void dust_list::setMaxYLoc(int maxY)
{
	maxYLoc = maxY;
}

int dust_list::newUniqueID()
{
	int currentID = uniqueID;
	uniqueID++;
	//std::cout << "Creating particle id: " << currentID <<  "." << std::endl;
	return currentID;
}

int dust_list::getVecLocByID(int id)
{
	for(int i = 0; i < myTotal; i++)
	{
		if( myDustList[i].getID() == id)
			return i;
	}
	std::cout << "Error locating dust grain " << id << "'s vector location." << std::endl;
	return -1;
}

dust_grain dust_list::getGrainByID(int id)   // gets element "id" from vector "myDustList"
{
	for (int i = 0; i < myTotal; i++)
	{
		if(myDustList[i].getID() == id)
			return myDustList[i];
	}
	std::cout << "Error copying dust grain, id " << id << ". Returning empty dust grain." << std::endl;
	return dust_grain();
}

void dust_list::incrementHandled()
{
	ptclsHandled++;
}

void dust_list::resetHandled()
{
	ptclsHandled = 0;
}

int dust_list::getIDByVecLoc(int n)
{
	return myDustList[n].getID();
}

void dust_list::dust_dstr()
{
	int s;
	std::vector <int> temp;
	bool add, filter;
	if (dustDist.size() == 0)
	{
		temp.push_back(myDustList[myTotal-1].getSize());
		temp.push_back(0);
		dustDist.push_back(temp);
		temp.clear();
	}
	for(unsigned int i = 0; i < myDustList.size(); i++)
	{
		add = false;
		s = myDustList[i].getSize();
		filter = myDustList[i].getFilter();
		if(!filter)
		{
			for(unsigned int j = 0; j < dustDist.size(); j++)
			{
				if(dustDist[j][0] == s)
				{
					dustDist[j][1]++;
					add = false;
					break;
				}
				else
					add = true;
			}
			if(add == true)
			{
				temp.push_back(s);
				temp.push_back(1);
				dustDist.push_back(temp);
				temp.clear();
			}
		}
	}
	temp.clear();
}

dust_grain dust_list::getGrainByVecLoc(int n)
{
	return myDustList[n];
}

void dust_list::update_width_dstr(int dWidth)
{
	bool add = false;
	std::vector <int> temp;

	if (dustWidth.size() == 0)
		add = true;
	else
	{
		for (unsigned int j = 0; j < dustWidth.size(); j++)
		{
			if (dustWidth[j][0] == dWidth)
			{
				dustWidth[j][1]++;
				add = false;
				break;
			}
			else
				add = true;
		}
	}
	if (add == true)
	{
		temp.push_back(dWidth);
		temp.push_back(1);
		dustWidth.push_back(temp);
		temp.clear();
	}
	temp.clear();
}

void dust_list::clear_width_dstr()
{
	dustWidth.clear();
}

std::vector <double> dust_list::calc_stats(std::vector< std::vector <int> > distribution)
{
	double pSum, sum, count, avg, std_dev, diff;
	count = 0;
	sum = 0;
	for (unsigned int i = 0; i < distribution.size(); i++)
	{
		pSum = distribution[i][0] * distribution[i][1];
		count = count + distribution[i][1];
		sum = sum + pSum;
	}
	avg = sum / count;

	for (unsigned int i = 0; i < distribution.size(); i++)
	{
		diff = avg - distribution[i][0];
		pSum = (diff * diff) * distribution[i][1];
		sum = sum + pSum;
	}

	std_dev = sum / count;
	std_dev = std::sqrt(std_dev);
	std::vector <double> stats(0);
	stats.push_back(avg);
	stats.push_back(std_dev);

	return stats;
}

void dust_list::update_dstr_merge(int old1, int old2, int n)
{
	bool add = false;
	std::vector <int> temp;
	for(unsigned int j = 0; j < dustDist.size(); j++)
	{
		if(dustDist[j][0] == n)
		{
			dustDist[j][1]++;
			add = false;
			break;
		}
		else
			add = true;
		if(dustDist[j][0] == old1 || dustDist[j][0] == old2)
		{
			if(dustDist[j][1] > 0)
				dustDist[j][1]--;
			else
				dustDist[j][1] = 0;
		}
	}
	if(add == true)
	{
		temp.push_back(n);
		temp.push_back(1);
		dustDist.push_back(temp);
		temp.clear();
	}
	temp.clear();
}

void dust_list::update_dstr_split(int old, int new1, int new2)
{
	bool add1 = false;
	bool add2 = false;
	bool done1 = false;
	bool done2 = false;
	std::vector <int> temp;
	for (unsigned int j = 0; j < dustDist.size(); j++)
	{
		if (dustDist[j][0] == new1 && !done1)
		{
			dustDist[j][1]++;
			add1 = false;
			done1 = true;
		}
		else if (!done1)
			add1 = true;

		if (dustDist[j][0] == new2 && !done2)
		{
			dustDist[j][1]++;
			add2 = false;
			done2 = true;
		}
		else if (!done2)
			add2 = true;

		if (dustDist[j][0] == old)
		{
			if (dustDist[j][1] > 0)
				dustDist[j][1]--;
			else
				dustDist[j][1] = 0;
		}
	}
	temp.clear();
	if (add1 == true)
	{
		temp.push_back(new1);
		temp.push_back(1);
		dustDist.push_back(temp);
		temp.clear();
	}
	if (add2 == true)
	{
		temp.push_back(new2);
		temp.push_back(1);
		dustDist.push_back(temp);
		temp.clear();
	}

}

void dust_list::update_dstr_stuck(int remove)
{
	for (unsigned int j = 0; j < dustDist.size(); j++)
	{
		if (dustDist[j][0] == remove)
		{
			if (dustDist[j][1] > 0)
				dustDist[j][1]--;
			else
				dustDist[j][1] = 0;
		}
	}
}

void dust_list::update_size_dstr(int dSize)
{
	bool add = false;
	std::vector <int> temp;

	if (sizeDist.size() == 0)
		add = true;
	else
	{
		for (unsigned int j = 0; j < sizeDist.size(); j++)
		{
			if (sizeDist[j][0] == dSize)
			{
				sizeDist[j][1]++;
				add = false;
				break;
			}
			else
				add = true;
		}
	}
	if (add == true)
	{
		temp.push_back(dSize);
		temp.push_back(1);
		sizeDist.push_back(temp);
		temp.clear();
	}
	temp.clear();
}

void dust_list::clear_size_dstr()
{
	sizeDist.clear();
}
