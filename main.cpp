//  Adam Sokolow
//  adam.sokolow@duke.edu
//  Dust simulation for Dr. Sen
//  Dated July 29 2005

/*	Edited by Kevin VanSlyke
kgvansly@buffalo.edu
Dated Jan 2 2016*/

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <iomanip>
#include <sstream>
#include <vector>
#include <math.h>
#include <errno.h>
#include <time.h>

/*	direct.h for windows filesystem manip
	unistd.h and sys/stat.h for unix filesystem manip*/
//#include <direct.h>
#include <unistd.h>
#include <sys/stat.h>
#include "timer.h"


#include "world.h"
#include "parameterReader.h"

using namespace std;


// *******
// M A I N
// *******
int main(int argc, char *argv[])
{
	int trialID;
	struct timespec before, after;
	struct timespec time_diff;
	double time_s;
	world * myWorld;
	int timeCount = 0;
	get_time(&before);


	/////////////////////
	/* Parameters so far:

	World Dimensions: X, Y
	World Speeds: x, y

	Total Dust: n
	Size Range of Dust: 1 -> N

	Total timesteps: t
	Enable: merging, splitting, sticking
	Filter Parameters: Length, Width of Slits, Width of Gaps
	trialID
	*/
	int filter, xSites, ySites, xMom, yMom, negYMom, totalGrains, minGrainSize, maxGrainSize, maxTime, sticking, splitting, merging, filterWidth, filterGap, filterLength, filter2Width, filter2Gap, filter2Length;
	bool enableSticking, enableSplitting, enableMerging;
	if (argc != 14 || argc != 17 || argc != 20)
	{
		std::cout << "Not enough parameters passed to main (Need 13, 16 or 19 integers after the executable name), trying to read from parameters.txt... " << std::endl;
		parameterReader *pR = new parameterReader();
		xSites = pR->getXSites();
		ySites = pR->getYSites();
		xMom = pR->getXMom();
		yMom = pR->getYMom();
		negYMom = pR->getNegYMom();
		totalGrains = pR->gettotalGrains();
		maxTime = pR->getMaxTime();
		enableSticking = pR->getSticking();
		enableSplitting = pR->getSplitting();
		enableMerging = pR->getMerging();
		minGrainSize = pR->getMinGrainSize();
		maxGrainSize = pR->getMaxGrainSize();
		filterWidth = pR->getFilterWidth();
		filterGap = pR->getFilterGap();
		filterLength = pR->getFilterLength();
		filter2Width = pR->getFilter2Width();
		filter2Gap = pR->getFilter2Gap();
		filter2Length = pR->getFilter2Length();
		trialID = pR->getTrialID();
	}
	else if(argc == 14)
	{
		xSites = atof(argv[1]);
		ySites = atof(argv[2]);
		xMom = atof(argv[3]);
		yMom = atof(argv[4]);
		negYMom = atof(argv[5]);
		totalGrains = atof(argv[6]);
		minGrainSize = atof(argv[7]);
		maxGrainSize = atof(argv[8]);
		maxTime = atof(argv[9]);
		sticking = atof(argv[10]);
		if (sticking == 1)
			enableSticking = true;
		else
			enableSticking = false;
		splitting = atof(argv[11]);
		if (splitting == 1)
			enableSplitting = true;
		else
			enableSplitting = false;
		merging = atof(argv[12]);
		if (merging == 1)
			enableMerging = true;
		else
			enableMerging = false;
		trialID = atof(argv[13]);
		filterWidth = -1;
		filterGap = -1;
		filterLength = -1;
		filter2Width = -1;
		filter2Gap = -1;
		filter2Length = -1;
		filter = 0;
	}
	else if(argc == 17)
	{
		xSites = atof(argv[1]);
		ySites = atof(argv[2]);
		xMom = atof(argv[3]);
		yMom = atof(argv[4]);
		negYMom = atof(argv[5]);
		totalGrains = atof(argv[6]);
		minGrainSize = atof(argv[7]);
		maxGrainSize = atof(argv[8]);
		maxTime = atof(argv[9]);
		sticking = atof(argv[10]);
		if (sticking == 1)
			enableSticking = true;
		else
			enableSticking = false;
		splitting = atof(argv[11]);
		if (splitting == 1)
			enableSplitting = true;
		else
			enableSplitting = false;
		merging = atof(argv[12]);
		if (merging == 1)
			enableMerging = true;
		else
			enableMerging = false;
		filterWidth = atof(argv[13]);
		filterGap = atof(argv[14]);
		filterLength = atof(argv[15]);
		trialID = atof(argv[16]);
		filter2Width = -1;
		filter2Gap = -1;
		filter2Length = -1;
		filter = 1;
	}
	else if(argc == 20)
	{
		xSites = atof(argv[1]);
		ySites = atof(argv[2]);
		xMom = atof(argv[3]);
		yMom = atof(argv[4]);
		negYMom = atof(argv[5]);
		totalGrains = atof(argv[6]);
		minGrainSize = atof(argv[7]);
		maxGrainSize = atof(argv[8]);
		maxTime = atof(argv[9]);
		sticking = atof(argv[10]);
		if (sticking == 1)
			enableSticking = true;
		else
			enableSticking = false;
		splitting = atof(argv[11]);
		if (splitting == 1)
			enableSplitting = true;
		else
			enableSplitting = false;
		merging = atof(argv[12]);
		if (merging == 1)
			enableMerging = true;
		else
			enableMerging = false;
		filterWidth = atof(argv[13]);
		filterGap = atof(argv[14]);
		filterLength = atof(argv[15]);
		filter2Width = atof(argv[16]);
		filter2Gap = atof(argv[17]);
		filter2Length = atof(argv[18]);
		trialID = atof(argv[19]);
		filter = 2;
	}

	myWorld = new world(xSites, ySites, xMom, yMom, negYMom);
	/*Sets the simulation size in the dust_list object for use in dust_list routines*/
	//TODO: Make arbitrary for any number of filter lines appended to end of parameter text file.
	if (filterGap == (-1) && filter2Gap == (-1))
	{
		filter = 0;
		myWorld->populateWorld(totalGrains, minGrainSize, maxGrainSize);
		std::cout << "No Filter detected." << std::endl;
	}
	else if (filterGap != (-1) && filter2Gap == (-1))
	{
		bool bimodal = false;
		filter = 1;
		if(bimodal)
		{
			myWorld->bimodalPopulateWorld(totalGrains, minGrainSize, maxGrainSize/4, 3*maxGrainSize/4, maxGrainSize, filterGap, filterWidth, filterLength);
		}
		else
			myWorld->populateWorld(totalGrains, minGrainSize, maxGrainSize, filterGap, filterWidth, filterLength);
		std::cout << "Filter gap: " << filterGap << std::endl;
	}
	else if (filterGap != (-1) && filter2Gap != (-1))
	{
		filter = 2;
		myWorld->populateWorld(totalGrains, minGrainSize, maxGrainSize, filterGap, filterWidth, filterLength, filter2Gap, filter2Width, filter2Length);
		std::cout << "Filter 1 gap: " << filterGap << std::endl;
		std::cout << "Filter 2 gap: " << filter2Gap << std::endl;
	}

	/*  Routine for creating folder, need to make alternate version for Linux file systems*/

	std::cout << "Total No. of Dust Grains =  " << myWorld->myList->getTotal() - filter << std::endl;
	std::cout << "X Sites = " << myWorld->getMaxXSize() << ". Y Sites = " << myWorld->getMaxYSize() << ". " << std::endl;
	//"/gpfs/scratch/kgvansly/"
	//"/projects/academic/sen/kgvansly/Dust_Data/"
	//"/home/kevin/Dust_Data/"
	std::ostringstream oFolder;
	oFolder << "/home/kevin/Dust_Data/" << filter << "fltrs" << filterGap << "pr" << filterWidth << "fbr" << filterLength << "fl"<< totalGrains << "ptcls" << minGrainSize << "-" << maxGrainSize << "dstr" << xSites << "x" << ySites << "y" << xMom << "Px" << yMom << "Py" << negYMom << "nPy" << maxTime << "tm";
	std::string outputFolder = oFolder.str();

	if(mkdir(outputFolder.c_str(), S_IRWXU) == -1)
		std::cout << "Folder " << outputFolder << " already exists, entering..." << std::endl;
	else
		mkdir(outputFolder.c_str(), S_IRWXU);

	struct stat fileInfo;
	std::string paramFile = outputFolder + "/parameters.txt";
	if (stat(paramFile.c_str(), &fileInfo) == 0)
	{
		std::cout << "parameters.txt already exists, not writting" << std::endl;
	}
	else
	{
		FILE * parameterFile = fopen(paramFile.c_str(), "a");
		// OUTPUT: dust, size, y-position, y-step, x-localSp, y-localSp
		fprintf(parameterFile, "%d %d \n%d %d %d\n%d \n%d %d \n%d \n%d %d %d \n%d %d %d \n%d %d %d", xSites, ySites, xMom, yMom, negYMom, totalGrains, minGrainSize, maxGrainSize, maxTime, sticking, splitting, merging, filterWidth, filterGap, filterLength, filter2Width, filter2Gap, filter2Length);
		fclose(parameterFile);
	}

	std::ostringstream pFolder;
	pFolder << outputFolder +  "/Trial" << trialID;
	std::string processFolder = pFolder.str();
	mkdir(processFolder.c_str(), S_IRWXU);
	myWorld->setProcOutputFolder(processFolder);
	myWorld->myList->setProcOutputFolder(processFolder);
	/* End of folder creation routine */
	/* Statistics collection calls, might be better to move them to lower level objects. */
	myWorld->myList->dust_dstr();
	myWorld->myList->setFunctionality(enableSplitting, enableSticking, enableMerging);
	//Debug line
	//std::cout << "Dust distribution has been calculated." << std::endl;

	/* Runs program for determined amount of time. */
	while (timeCount < maxTime)
	{
		if (timeCount % 500 == 0)
			std::cout << "Time in world  = " << timeCount << " " << std::endl;
		timeCount = timeCount + 1;
		myWorld->takeStep();
		myWorld->updateWorld();
	}
	get_time(&after);
	diff(&before, &after, &time_diff);
	time_s = time_diff.tv_sec + (double)(time_diff.tv_nsec)/1.0e9;
	std::cout << "Total runtime is: " << time_s << " seconds." << std::endl;
	std::string timeOptFile = outputFolder + "/timeOptimization.txt";
	ofstream timeOptimization(timeOptFile.c_str(), ios::app);
	if (timeOptimization.is_open())
	{
		//writing: time of calculation, number of particles and length/width
		timeOptimization << "Total runtime is: " << time_s << " seconds." << std::endl;
		timeOptimization << "For " << maxTime << "timesteps." << std::endl;
		timeOptimization << "Simulation of " << totalGrains << " particles in an area of:" << std::endl;
		timeOptimization << xSites << " by " << ySites << ". (x by y)" << std::endl;
		timeOptimization.close();
	}

	std::cout << "Simulation completed!" << std::endl;
	return 0;
}
