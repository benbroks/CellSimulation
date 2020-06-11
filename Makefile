all: simulation

simulation: simulation.cpp cell.cpp colony.cpp
	g++ -g -Wall simulation.cpp cell.cpp colony.cpp -o simulation

clean:
	rm simulation
	rm -rf simulation.dSYM
	rm .DS_Store