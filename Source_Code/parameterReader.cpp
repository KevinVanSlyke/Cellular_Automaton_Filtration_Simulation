
//  Adam Sokolow
//  adam.sokolow@duke.edu
//  Dust simulation for Dr. Sen
//  Dated July 29 2005

/*	Edited by Kevin VanSlyke
kgvansly@buffalo.edu
Dated Feb 7 2016	*/

#include <ctime>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include "parameterReader.h"

//Loads in manditory parameters, with options that allow for loading parameters for up to two distict filters.
parameterReader::parameterReader()
{
	std::ifstream in_file("parameters.txt");
	std::string line;
	int lnNum = 1;
	while(std::getline(in_file, line))
	{
		std::istringstream iss(line);
		switch (lnNum)
		{
			case 1:
			{
				iss >> xSites >> ySites;
				std::cout << "Line " << lnNum << " reads: xSites = " << xSites << ", ySites = " << ySites << "." << std::endl;
				break;
			}
			case 2:
			{
				iss >> xMom >> yMom >> negYMom;
				std::cout << "Line " << lnNum << " reads: xMom = " << xMom << ", yMom = " << yMom << ", negYMom = " << negYMom << "." << std::endl;
				break;
			}
			case 3:
			{
				iss >> totalGrains;
				std::cout << "Line " << lnNum << " reads: totalGrains = " << totalGrains << "." << std::endl;
				break;
			}
			case 4:
			{
				iss >> minGrainSize >> maxGrainSize;
				std::cout << "Line " << lnNum << " reads: minGrainSize = " << minGrainSize << ", maxGrainSize = " << maxGrainSize << "." << std::endl;
				break;
			}
			case 5:
			{
				iss >> maxTime;
				std::cout << "Line " << lnNum << " reads: maxTime = " << maxTime << "." << std::endl;
				break;
			}
			case 6:
			{
				iss >> sticking >> splitting >> merging;
				std::cout << "Line " << lnNum << " reads: sticking = " << sticking << ", splitting = " << splitting << ", merging = " << merging << "." << std::endl;
				break;
			}
			case 7:
			{
				iss >> filterWidth >> filterGap >> filterLength;
				std::cout << "Line " << lnNum << " reads: filterWidth = " << filterWidth << ", filterGap = " << filterGap << ", filterLength = " << filterLength << "." << std::endl;
				if(filterWidth < 0 || filterGap < 0 || filterLength < 0)
					filter = false;
				else
					filter = true;
				break;
			}
			case 8:
			{
				iss >> trialID;
				std::cout << "Line " << lnNum << " reads: trialID = " << trialID << "." << std::endl;
				break;
			}
			case 9:
			{
				iss >> filter2Width >> filter2Gap >> filter2Length;
				std::cout << "Line " << lnNum << " reads: filter2Width = " << filter2Width << ", filter2Gap = " << filter2Gap << ", filter2Length = " << filter2Length << "." << std::endl;
				filter2 = true;
				break;
			}

		}
	++lnNum;
	}
	in_file.close();
}

int parameterReader::getMaxTime()
{
	return maxTime;
}

bool parameterReader::getSticking()
{
	if (sticking == 1)
		return true;
	else
		return false;
}

bool parameterReader::getSplitting()
{
	if (splitting == 1)
		return true;
	else
		return false;
}

bool parameterReader::getMerging()
{
	if (merging == 1)
		return true;
	else
		return false;
}

int parameterReader::getXSites()
{
	return xSites;
}

int parameterReader::getYSites()
{
	return ySites;
}

int parameterReader::getXMom()
{
	return xMom;
}
int parameterReader::getYMom()
{
	return yMom;
}
int parameterReader::getNegYMom()
{
	return negYMom;
}
int parameterReader::gettotalGrains() //Credited Adam Sokolow
{
	return totalGrains;
}

int parameterReader::getMaxGrainSize() //Credited Adam Sokolow
{
	return maxGrainSize;
}

int parameterReader::getMinGrainSize() //Credited Adam Sokolow
{
	return minGrainSize;
}

int parameterReader::getFilterWidth() //Credited Adam Sokolow
{
	if (filter)
		return filterWidth;
	else
		return -1;
}

int parameterReader::getFilterGap() //Credited Adam Sokolow
{

	if (filter)
		return filterGap;
	else
		return -1;
}
int parameterReader::getFilterLength() //Credited Adam Sokolow
{

	if (filter)
		return filterLength;
	else
		return -1;
}

int parameterReader::getFilter2Width()
{
	if (filter2)
		return filter2Width;
	else
		return -1;
}

int parameterReader::getFilter2Gap()
{

	if (filter2)
		return filter2Gap;
	else
		return -1;
}
int parameterReader::getFilter2Length()
{

	if (filter2)
		return filter2Length;
	else
		return -1;
}
int parameterReader::getTrialID()
{
	return trialID;
}