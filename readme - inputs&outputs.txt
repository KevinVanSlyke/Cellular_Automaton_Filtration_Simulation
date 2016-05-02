Inputs:
i) parameters.txt
Parameters are read from top to bottom as follows:
xMax yMax (horizontal/width of simulation area) (vertical/height of simulation area)
xSpeed ySpeed (maximum movement of particle in the given dimension)
totalGrains 
MinGrainSize MaxGrainSize
maxTimeSteps
enable: sticking, splitting, merging
trialID
FilterWidth FilterGap FilterLength
Filter2Width Filter2Gap Filter2Length

Outputs:
1) dustCount.txt
At every timestep: # Dust moving, stuck, too large to move, merged, and total

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



