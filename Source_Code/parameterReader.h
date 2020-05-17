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
	int getXSites();
	int getYSites();
	int getXMom();
	int getYMom();
	int getNegYMom();
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
	int xSites;
	int ySites;
	int xMom;
	int yMom;
	int negYMom;
	int totalGrains;
	int maxGrainSize;
	int minGrainSize;
	int maxTime;
	int filterWidth;
	int filterGap;
	int filterLength;
	int filter2Width;
	int filter2Gap;
	int filter2Length;
	bool filter;
	bool filter2;
	int sticking;
	int splitting;
	int merging;
	int trialID;
};
#define parameterReader_H 
#endif 
