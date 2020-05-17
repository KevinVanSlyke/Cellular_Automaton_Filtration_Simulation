
//  Adam Sokolow 
//  adam.sokolow@duke.edu 
//  Dust simulation for Dr. Sen 
//  Dated July 29 2005 

/*	Edited by Kevin VanSlyke
kgvansly@buffalo.edu
Dated Jan 2 2016 */

#ifndef DUST_GRAIN_H 



#include<vector> 
#include<iostream> 
#include<cstdlib> 



// ******************************************************************* 
// Create a dust_grain Class that keeps track of the visual location of the dust 

// ******************************************************************* 

class dust_grain
{
public:

	// Constructors/Destructor 
	dust_grain();                        // default constructor (size==500)
	dust_grain(std::vector<int> x, std::vector<int> y, int size);
	dust_grain(std::vector<int> x, std::vector<int> y, int size, int id);    // initial size of dust_grain 
	dust_grain(const dust_grain & d);   // copy constructor 
	~dust_grain();                       // destructor 

  	// Assignment 
	const dust_grain & operator = (const dust_grain & d);

	// Accessors
	bool checkMoved();
	bool checkMerge();
	bool checkSplit();
	bool getFilter();
	bool getStuck();
	bool spotTaken(int x, int y);
	int getSize();
	int getXatc(int c);
	int getYatc(int c);
	int getPrevYMom();
	int getPrevXMom();
	int getID();
	int getPrevPB();
	int getCurPB();
	int getMaxXStep();
	int getMaxYStep();
	int getColXMom();
	int getColYMom();
	std::vector <int> getXent();
	std::vector <int> getYent();
	
	// Simple modifiers
	void setSize();
	void setPrevYMom(int yMom);
	void setPrevXMom(int xMom);
	void setID(int newID);
	void setStuck(bool stk);
	void setPrevPB(int num);
	void setCurPB(int num);
	void setFilter(bool filt);
	void setMoved(bool moved);
	void setMerge(bool merge);
	void setSplit(bool split);
	void setMaxXLoc(int xLen);
	void setMaxYLoc(int yLen);
	void setMaxXStep(int mXStep);
	void setMaxYStep(int mYStep);
	void setColXMom(int cXMom);
	void setColYMom(int cYMom);
	// Complex Modifiers
	void moveStep(int x, int y);
	void growGrain(std::vector <int> x, std::vector <int> y);
	int calculateWidth();
	
	// Clear data
	void clearGrain();
	
	//std::vector < dust_grain > attemptBreakUp(int maxX, int maxY)

private:
	bool stuck;
	bool pendingMerge;
	bool pendingSplit;
	bool hasMoved;
	bool filter;
	
	int mySize;  //dust_grain Size 
	int prevPillbox;
	int curPillbox;
	int grainID;
	int maxXLoc;
	int maxYLoc;
	int width;
	int maxXStep;
	int maxYStep;
	int prevXMom;
	int prevYMom;
	int colXMom;
	int colYMom;
	std::vector<int> myX, myY;
};
#define DUST_GRAIN_H 
#endif 
