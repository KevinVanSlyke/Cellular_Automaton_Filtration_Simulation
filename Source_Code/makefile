# the compiler: gcc for C program, define as g++ for C++
CC = g++
CFLAGS = -g -Wall -O0
LIBS = -lGL -lGLU -lglut
# the build target executable:
all: dust_simulation visual_dust_simulation

visual_dust_simulation : visualMain.cpp world.cpp dustlist.cpp dustgrain.cpp parameterReader.cpp randomgen.cpp timer.cpp
	$(CC) $(CFLAGS) -o visual_dust_simulation visualMain.cpp world.cpp dustlist.cpp dustgrain.cpp parameterReader.cpp randomgen.cpp timer.cpp $(LIBS)

dust_simulation : main.cpp world.cpp dustlist.cpp dustgrain.cpp parameterReader.cpp randomgen.cpp timer.cpp
	$(CC) $(CFLAGS) -o dust_simulation main.cpp world.cpp dustlist.cpp dustgrain.cpp parameterReader.cpp randomgen.cpp timer.cpp

clean:
	-rm -f ./dust_simulation ./visual_dust_simulation
