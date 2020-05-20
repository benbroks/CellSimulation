all: simulation

simulation: simulation.cpp cell.cpp
	g++ -g -Wall simulation.cpp cell.cpp -o simulation

clean:
	rm simulation
	rm -rf simulation.dSYM