
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

	// constructors/destructor 
	dust_grain();                        // default constructor (size==500)
	dust_grain(std::vector<int> x, std::vector<int> y, int size);
	dust_grain(std::vector<int> x, std::vector<int> y, int size, int id);    // initial size of dust_grain 
	dust_grain(const dust_grain & d);   // copy constructor 
	~dust_grain();                       // destructor 

  // assignment 
	const dust_grain & operator = (const dust_grain & d);

	// accessors 
	int getSize();
	void setSize();
	bool spotTaken(int x, int y);
	int getXatc(int c);
	int getYatc(int c);
	std::vector <int> getXent();
	std::vector <int> getYent();
	// modifiers

	//TODO: Should make each grain initialize with the value of maxX/Y Loc or figure out how to have it access it's parent.
	void moveStep(int x, int y, int maxX, int maxY);
	void setStuck(bool stk);
	bool getStuck();
	int getPrevPB();
	int getCurPB();
	void setPrevPB(int num);
	void setCurPB(int num);
	void setFilter(bool filt);
	bool getFilter();
	
	int getID();
	void setID(int newID);
	bool checkMoved();
	bool checkMerge();
	void setMoved(bool moved);
	void setMerge(bool merge);
	void setMaxXVel(int xvel);
	int getMaxXVel();
	void setMaxYVel(int yvel);
	int getMaxYVel();
	void growGrain(std::vector <int> x, std::vector <int> y);
	void clearGrain();
	int calculateWidth(int maxXLoc);
	//std::vector < dust_grain > attemptBreakUp(int maxX, int maxY)

private:

	int mySize;  //dust_grain Size 
	std::vector<int> myX, myY;
	bool stuck;
	int prevPillbox;
	int curPillbox;
	
	int grainID;
	bool pendingMerge;
	bool hasMoved;
	int maxXVel;
	int maxYVel;
	bool filter;
	int width;
};
#define DUST_GRAIN_H 
#endif 
