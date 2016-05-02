
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
				iss >> xMax >> yMax;
				std::cout << "Line " << lnNum << " reads: xMax = " << xMax << ", yMax = " << yMax << "." << std::endl;
				break;
			}
			case 2:
			{
				iss >> xSpeed >> ySpeed;
				std::cout << "Line " << lnNum << " reads: xSpeed = " << xSpeed << ", ySpeed = " << ySpeed << "." << std::endl;
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
				iss >> MinGrainSize >> MaxGrainSize;
				std::cout << "Line " << lnNum << " reads: MinGrainSize = " << MinGrainSize << ", MaxGrainSize = " << MaxGrainSize << "." << std::endl;
				break;
			}
			case 5:
			{
				iss >> MaxTime;
				std::cout << "Line " << lnNum << " reads: MaxTime = " << MaxTime << "." << std::endl;
				break;
			}
			case 6:
			{
				iss >> sticking >> splitting >> merging;
				std::cout << "Line " << lnNum << " reads: Sticking = " << sticking << ", Splitting = " << splitting << ", Merging = " << merging << "." << std::endl;
				break;
			}
			case 7:
			{
				iss >> trialID;
				std::cout << "Line " << lnNum << " reads: TrialID = " << trialID << "." << std::endl;
				break;
			}
			case 8:
			{
				iss >> FilterWidth >> FilterGap >> FilterLength;
				std::cout << "Line " << lnNum << " reads: FilterWidth = " << FilterWidth << ", FilterGap = " << FilterGap << ", FilterLength = " << FilterLength << "." << std::endl;
				Filter = true;
				break;
			}
			case 9:
			{
				iss >> Filter2Width >> Filter2Gap >> Filter2Length;
				std::cout << "Line " << lnNum << " reads: Filter2Width = " << Filter2Width << ", Filter2Gap = " << Filter2Gap << ", Filter2Length = " << Filter2Length << "." << std::endl;
				Filter2 = true;
				break;
			}

		}
	++lnNum;
	}
	in_file.close();
}

int parameterReader::getMaxTime()
{
	return MaxTime;
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

int parameterReader::getxMax()
{
	return xMax;
}

int parameterReader::getyMax()
{
	return yMax;
}

int parameterReader::getxSpeed()
{
	return xSpeed;
}
int parameterReader::getySpeed()
{
	return ySpeed;
}
int parameterReader::gettotalGrains()
{
	return totalGrains;
}

int parameterReader::getMaxGrainSize()
{
	return MaxGrainSize;
}

int parameterReader::getMinGrainSize()
{
	return MinGrainSize;
}

int parameterReader::getFilterWidth()
{
	if (Filter)
		return FilterWidth;
	else
		return -1;
}

int parameterReader::getFilterGap()
{

	if (Filter)
		return FilterGap;
	else
		return -1;
}
int parameterReader::getFilterLength()
{

	if (Filter)
		return FilterLength;
	else
		return -1;
}

int parameterReader::getFilter2Width()
{
	if (Filter2)
		return Filter2Width;
	else
		return -1;
}

int parameterReader::getFilter2Gap()
{

	if (Filter2)
		return Filter2Gap;
	else
		return -1;
}
int parameterReader::getFilter2Length()
{

	if (Filter2)
		return Filter2Length;
	else
		return -1;
}
int parameterReader::getTrialID()
{
	return trialID;
}