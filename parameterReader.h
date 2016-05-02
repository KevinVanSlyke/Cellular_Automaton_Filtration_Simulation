//  Adam Sokolow 
//  adam.sokolow@duke.edu 
//  Dust simulation for Dr. Sen 
//  Dated July 29 2005 

/*	Edited by Kevin VanSlyke
kgvansly@buffalo.edu
Dated Jan 2 2016 */

#ifndef PARAMETERREADER_H 

#include<iostream> 
#include<fstream> 
#include<cstdlib>
#include<vector>

class parameterReader
{
public:
	parameterReader();
	int getxMax();
	int getyMax();
	int getxSpeed();
	int getySpeed();
	int gettotalGrains();
	int getMaxGrainSize();
	int getMinGrainSize();
	int getMaxTime();
	bool getSticking();
	bool getSplitting();
	bool getMerging();
	int getFilterWidth();
	int getFilterGap();
	int getFilterLength();
	int getFilter2Width();
	int getFilter2Gap();
	int getFilter2Length();
	int getTrialID();

private:
	int xMax;
	int yMax;
	int xSpeed;
	int ySpeed;
	int totalGrains;
	int MaxGrainSize;
	int MinGrainSize;
	int MaxTime;
	int FilterWidth;
	int FilterGap;
	int FilterLength;
	int Filter2Width;
	int Filter2Gap;
	int Filter2Length;
	bool Filter;
	bool Filter2;
	int sticking;
	int splitting;
	int merging;
	int trialID;
};
#define parameterReader_H 
#endif 
