./Source_Code contains the .cpp, .h and makefile to compile the agent based particle filtration software.
./Executable contains stable executable with and without visualization, along with scripts to run them.
	Note that every time the makefile is run in ./Source_Code the new exectuables need to be moved to ./Executable
./Executable/parameters.txt is read by the executable in the following format to set parameters.

Inputs:
i) parameters.txt
Parameters are read from top to bottom as follows:
xSites ySites (horizontal/width of simulation area in # of lattice sites) (vertical/height of simulation area)
xMom yMom yNegMom (max momomentum of particle in the given dimension that can be randomly caused by the fluid, min x mom = -xMom, third option tells reverse y momentum)
totalGrains
MinGrainSize MaxGrainSize
maxTimeSteps
enable: sticking, splitting, merging
trialID
poreSeparation poreWidth poreDepth (for filter layer 1)
pore2Separation pore2Width pore2Depth (for filter layer 2)

Outputs:
1) dustCount.txt
At every timestep: # Dust moving, stuck, pseudo-stuck (not moving near filter), too large to move, merged, total, handled

2) dustDist.txt
Distribution of sizes and widths
time, avg width of moving ptcls, width std dev of moving ptcls, avg size of moving ptcls, size std dev of moving ptcls, avg size of all ptcls, size std dev of all ptcls

3) dustfile.txt
This file is huge. At every time:
grainID, size, x-position of pixel-0, y-position of pixel-0, curx, cury

4) dustfilePillCount.txt
At each timestep:
# of ptcls within pillBox around each pore

5) poresBlocked.txt
ID of pore blocked, time of blockage
