//  Adam Sokolow
//  adam.sokolow@duke.edu
//  Dust simulation for Dr. Sen
//  Dated July 29 2005

/*	Edited by Kevin VanSlyke
kgvansly@buffalo.edu
Dated Jan 2 2016 */

#include <ctime>
#include <cstdlib>
#include "dustgrain.h"
#include <malloc.h>
#include <vector>

//Default constructor
dust_grain::dust_grain()
{
	mySize = 0;
	myX = std::vector<int>(0);
	myY = std::vector<int>(0);
	stuck = false;
	prevPillbox = -1;
	curPillbox = -1;
	pendingMerge = false;
	hasMoved = false;
	grainID = -1;
	filter = false;
}

//Constructor with parameters
dust_grain::dust_grain(std::vector<int> x, std::vector<int> y, int size)
{
	mySize = size;
	myX = x;
	myY = y;
	stuck = false;
	prevPillbox = -1;
	curPillbox = -1;
	pendingMerge = false;
	hasMoved = false;
	grainID = -1;
	filter = false;
}

dust_grain::dust_grain(std::vector<int> x, std::vector<int> y, int size, int id)
{
	mySize = size;
	myX = x;
	myY = y;
	stuck = false;
	prevPillbox = -1;
	curPillbox = -1;
	pendingMerge = false;
	hasMoved = false;
	grainID = id;
	filter = false;
	//Debug Line
	//std::cout << "Grain ID " << grainID << " created." << std::endl;
}

//Copier.
dust_grain::dust_grain(const dust_grain  & d)
{

	mySize = d.mySize;
	myX = d.myX;
	myY = d.myY;
	stuck = d.stuck;
	prevPillbox = d.prevPillbox;
	curPillbox = d.curPillbox;
	pendingMerge = d.pendingMerge;
	hasMoved = d.hasMoved;
	grainID = d.grainID;
	filter = d.filter;
}

//Desctructor.
dust_grain ::~dust_grain()
{
	mySize = 0;
}

const dust_grain  &
dust_grain ::operator = (const dust_grain  & rhs)
{
	if (this != &rhs)                           // don't assign to self!
	{

		mySize = rhs.mySize;
		myX = rhs.myX;
		myY = rhs.myY;
		stuck = rhs.stuck;
		prevPillbox = rhs.prevPillbox;
		curPillbox = rhs.curPillbox;
		pendingMerge = rhs.pendingMerge;
		hasMoved = rhs.hasMoved;
		grainID = rhs.grainID;
		filter = rhs.filter;
	}
	return *this;
}

int dust_grain::getSize()
{
	return mySize;
}

void dust_grain::setSize()
{
	mySize = myX.size();
}

bool dust_grain::getStuck()
{
	return stuck;
}

void dust_grain::setStuck(bool stk)
{
	stuck = stk;
}

//Check if pixel is occupied by itself.
bool dust_grain::spotTaken(int x, int y)
{
	for (int c = 0; c < mySize; c++)
		if (myX[c] == x && myY[c] == y)
			return true;

	return false;
}

//Moves the location of a single dust particle.
void dust_grain::moveStep(int x, int y, int maxX, int maxY)
{
	for (int c = 0; c < mySize; c++)
	{
		myX[c] = (myX[c] + x + maxX) % maxX;
		myY[c] = (myY[c] + y + maxY) % maxY;
	}
}

int dust_grain::getXatc(int c)
{
	if (c < mySize)
		return myX[c];
	else
		return -1;
}

int dust_grain::getYatc(int c)
{

	if (c < mySize)
		return myY[c];
	else
		return -1;
}

int dust_grain::getPrevPB()
{
	return prevPillbox;
}

void dust_grain::setPrevPB(int num)
{
	prevPillbox = num;
}

int dust_grain::getCurPB()
{
	return curPillbox;
}

void dust_grain::setCurPB(int num)
{
	curPillbox = num;
}

int dust_grain::getID()
{
	return grainID;
}

bool dust_grain::checkMoved()
{
	return hasMoved;
}

bool dust_grain::checkMerge()
{
	return pendingMerge;
}

void dust_grain::setMoved(bool moved)
{
	hasMoved = moved;
}

void dust_grain::setMerge(bool merge)
{
	pendingMerge = merge;
}

void dust_grain::setMaxXVel(int xvel)
{
	maxXVel = xvel;
}

int dust_grain::getMaxXVel()
{
	return maxXVel;
}
void dust_grain::setMaxYVel(int yvel)
{
	maxYVel = yvel;
}

int dust_grain::getMaxYVel()
{
	return maxYVel;
}

void dust_grain::setFilter(bool filt)
{
	filter = filt;
}

bool dust_grain::getFilter()
{
	return filter;
}

void dust_grain::setID(int newID)
{
	grainID = newID;
}

void dust_grain::growGrain(std::vector <int> x, std::vector <int> y)
{
	myX = x;
	myY = y;
	mySize = x.size();
}

void dust_grain::clearGrain()
{
	myX.clear();
	myY.clear();
	mySize = 0;
}

std::vector <int> dust_grain::getXent()
{
	return myX;
}

std::vector <int> dust_grain::getYent()
{
	return myY;
}

int dust_grain::calculateWidth(int maxXLoc)
{
	int hX = 0;
	int lX = maxXLoc * 2;
	int xPos;
	std::vector < int > copyMyX = std::vector < int >(mySize);

	//	bool boundary;

	for (int i = 0; i < mySize; i++)
	{
		xPos = getXatc(i);
		copyMyX[i] = xPos;
		if (xPos < lX)
			lX = xPos;
		if (xPos > hX)
			hX = xPos;

		//Debug statement to make sure that particles at the boundaries are handled properly

//		if (xPos <= 1)
//			boundary = true;
//		else
//			boundary = false;

	}
//	if(boundary)
//		std::cout << "Particle " << getID() << " is lying on a boundary during width measurement." << std::endl;
	bool xEdge = false;
	int xDiff = hX - lX;

	if (xDiff >= mySize)
	{
		xEdge = true;
		for (int i = 0; i < mySize; i++)
		{
			if (copyMyX[i] < (mySize - 1))
			{
				copyMyX[i] = copyMyX[i] + maxXLoc;
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
		xDiff = hX - lX;
	}
	width = xDiff+1;
	return width;
}

/*	//Splits a grain vertically or horizontally at a central cell, whichever is more even
//Has some major issues which caused me to pursue a better method.
std::vector < dust_grain > dust_grain::attemptBreakUp(int maxX, int maxY)
{
	std::vector < int > copyMyX = std::vector < int >(mySize);
	std::vector < int > copyMyY = std::vector < int >(mySize);
	std::vector < int > xDiffList = std::vector < int >(mySize);
	std::vector < int > yDiffList = std::vector < int >(mySize);
	std::vector < int > xSameList = std::vector < int >(mySize);
	std::vector < int > ySameList = std::vector < int >(mySize);
	int tempLoc, tempX, tempY, tempDiff, counter;
	bool vertCut;

	//Create a copy of continuous dust occupied cells to account for a grain being partially across a boundary
	for (int i = 0; i < mySize; i++)
	{
		copyMyX.push_back(myX[i] + maxX);
		copyMyY.push_back(myY[i] + maxY);
	}
	counter = 0;
	while (counter < mySize)
	{
		tempX = copyMyX[counter];
		tempY = copyMyY[counter];
		for (int i = 0; i < mySize; i++)
		{
			if (copyMyX[i] < tempX)
			{
				xDiffList[counter]--;
			}
			if (copyMyY[i] < tempY)
			{
				yDiffList[counter]--;
			}
			if (copyMyX[i] > tempX)
			{
				xDiffList[counter]++;
			}
			if (copyMyY[i] > tempY)
			{
				yDiffList[counter]++;
			}
			if (copyMyX[i] == tempX)
			{
				xSameList[counter]++;
			}
			if (copyMyY[i] == tempY)
			{
				ySameList[counter]++;
			}
		}
		counter++;
	}
	tempDiff = std::abs(xDiffList[0]);
	for ( int i = 0; i < mySize; i++)
	{
		if ((xDiffList[i] == 0) && (xSameList[i] % 2 == 0))
		{
			tempLoc = i;
			vertCut = true;
			break;
		}
		else if ((yDiffList[i] == 0) && (ySameList[i] % 2 == 0))
		{
			tempLoc = i;
			vertCut = false;
			break;
		}
		else if (xDiffList[i] == 0)
		{
			tempLoc = i;
			tempDiff = std::abs(xDiffList[i]);
			vertCut = true;
		}
		else if (yDiffList[i] == 0)
		{
			tempLoc = i;
			tempDiff = std::abs(yDiffList[i]);
			vertCut = false;
		}
		else if (std::abs(xDiffList[i]) < tempDiff)
		{
			tempLoc = i;
			tempDiff = std::abs(xDiffList[i]);
			vertCut = true;
		}
		else if (std::abs(yDiffList[i]) < tempDiff)
		{
			tempLoc = i;
			tempDiff = std::abs(yDiffList[i]);
			vertCut = false;
		}
	}
	std::vector < int > firstShardX;
	std::vector < int > firstShardY;
	std::vector < int > secondShardX;
	std::vector < int > secondShardY;
	int gCount = 0;
	int lCount = 0;
	tempX = myX[tempLoc];
	tempY = myY[tempLoc];
	if (vertCut)
	{
		for (int i = 0; i < mySize; i++)
		{
			if (copyMyX[i] > tempX)
			{
				firstShardX.push_back(copyMyX[i]);
				firstShardY.push_back(copyMyY[i]);
			}
			else if (copyMyX[i] < tempX)
			{
				secondShardX.push_back(copyMyX[i]);
				secondShardY.push_back(copyMyY[i]);
			}
			else if (copyMyX[i] = tempX)
			{
				firstShardX.push_back(copyMyX[i]);
				firstShardY.push_back(copyMyY[i]);
			}
		}
	}
}	*/
