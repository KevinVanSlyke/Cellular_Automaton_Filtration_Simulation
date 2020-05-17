//  Adam Sokolow 
//  adam.sokolow@duke.edu 
//  Dust simulation for Dr. Sen 
//  Dated July 29 2005 

/*	Edited by Kevin VanSlyke
kgvansly@buffalo.edu
Dated Jan 2 2016 */

#ifndef DUST_LIST_H 

#include <vector> 
#include <iostream> 
#include <cstdlib> 
#include "randomgen.h" 
#include "dustgrain.h" 

// ******************************************************************* 
// Create a dust_list Class that keeps track of the visual location of the dust 

// ******************************************************************* 

class dust_list {
public:

	// constructors/destructor 
	dust_list();
	dust_list(int**   world);
	// default constructor (size==500) 
	dust_list(int**   world, int total, int low, int high);
	// initial size of dust_list 
	dust_list(const dust_list & dl);   // copy constructor 
	~dust_list();                       // destructor 

	// assignment 
	const dust_list & operator = (const dust_list & dl);

	// accessors 
	int getTotal();

	double Abs(double Nbr);
	int getTimeSteps();

	std::vector< std::vector < int > > getPillBoxes();
	std::vector < int > getPoreJamTimer();
	std::vector < bool > getPotentialBlock();
	int getpBCounts(int loc);
	bool getPoreBlocked(int loc);

	int getMaxXLoc();
	int getMaxYLoc();
	int getMaxXMom();
	int getMaxYMom();
	int getNegYMom();

	void setPillBoxes(int loc, std::vector< int > pBoxes);
	void setPBCounts(int loc, int pCounts);
	void setPoreJamTimer(int loc, int jamTimer);
	void setPoreBlocked(int loc, bool verBlock);
	void setPotentialBlock(int loc, bool potBlock);

	void setMaxXLoc(int maxX);
	void setMaxYLoc(int maxY);
	void setMaxXMom(int yMom);
	void setMaxYMom(int yMom);
	void setNegYMom(int negYMom);
	void resetPBCounts();
	void resetPotentialBlock();
	void incrimentPBCounts(int loc);
	void incrimentPoreJamCounter(int loc);
	void incrimentTimeStep();
	
	// modifiers 
	void moveStep(int** &updateWorld);

	void addGrain(int low, int high);
	void addGrain(int filterGap, int filterWidth, int filterLength);
	void addGrain2(int filter2Gap, int filter2Width, int filter2Length);
	//void addGrainSk(int filterGap, int filterWidth, int filterLength, float timeF, int ** &updateWorld);

	dust_grain getGrainByID(int id);
	dust_grain getGrainByVecLoc(int n);
	int getIDByVecLoc(int n);
	int newUniqueID();			//Returns an int then incriments uniqueID
	int getVecLocByID(int id);

	void dust_dstr();
	void setNewTotal();

	void setFunctionality(bool splitting, bool sticking, bool merging);

	void setProcOutputFolder(std::string dirName);

private:
	void removeMergedSplitGrains();
	void addMergedSplitGrains();

	void shrinkListbyOne();
	void increaseListbyOne();

	bool isOpen(int x, int y);
	bool isOpen(int x, int y, int ignore_grain_j);
	bool isOpenSelf(std::vector<int> x, std::vector<int> y, int ranSite, int ranCardinal, int psz);
	
	bool canMakeMove(int xmove, int ymove, int grainNumber);
	bool canMakeMove(int xmove, int ymove, dust_grain cgrain);

	int getCollidingGrain(int xmove, int ymove, int grain_self);
	//dust_grain mergeGrain_to_filter(int g1, int g2);
	dust_grain mergeGrains(std::vector<dust_grain> grainMergeSet);


	std::vector<dust_grain> splitGrain(dust_grain aGrain);
	void separateSplitGrains(dust_grain g1, dust_grain g2);

	std::vector<int> createRandomOrder();

	//Functions for creating pores and keeping track of if pores are blocked
	void calcPillBoxes(int filterGap, int filterWidth, int filterLength, int start);
	void checkBlocked();
	std::vector < bool > checkPoreFilled();
	void pushBackPillBoxVectors(std::vector < int > corners);

	//Older counting method, still used to track all particles
	void update_dstr_merge(std::vector<int> oldSizes, int n);
	void update_dstr_split(int old, int new1, int new2);
	void update_dstr_stuck(int remove);
	//For tracking moving particles area/size
	void update_size_dstr(int dSize);
	void clear_size_dstr();
	//For tracking moving particles widths
	void update_width_dstr(int dWidth);
	void clear_width_dstr();
	//Function to calculate size/width statistics
	std::vector <double> calc_stats(std::vector < std::vector < int > > distribution);

	void calculatePostCollisionMomentum(dust_grain movingPtcl, dust_grain staticPtcl, int movXMomBefore, int movYMomBefore);
	//Private variables
	int myTotal;
	int numTimeSteps;
	int maxXLoc;
	int maxYLoc;
	int maxXMom;
	int maxYMom;
	int negYMom;
	bool enableSticking;
	bool enableMerging;
	bool enableSplitting;
	random_gen* myGenerator;

	std::vector < dust_grain > myDustList;
	std::vector < dust_grain > grainsToAdd;
	std::vector < std::vector < dust_grain > > grainsToMerge;

	void addGrainsToMergeList(dust_grain g1, dust_grain g2);

	//begin new stuff to fix splitting and merging
/*	int pendingTotal;
	std::vector < dust_grain > pendingDustList;
	void syncPendingDustList();
	void updateMyDustList();

	int** pendingWorld;
	void syncPendingWorld();
	void updateRefWorld();

	//end new stuff
*/
	int** refWorld;

	int uniqueID;

	std::vector< std::vector< int > > pillBoxes;
	std::vector< int > pBCounts;
	std::vector< int > poreJamTimer;
	std::vector< bool > poreBlocked;
	std::vector< bool > potentialBlock;

	std::vector < std::vector < int > > dustDist;
	std::vector < std::vector < int > > sizeDist;
	std::vector < std::vector < int > > dustWidth;

	std::string procOutputFolder;

};
#define DUST_LIST_H 
#endif 

