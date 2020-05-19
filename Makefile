all: simulation

simulation: simulation.cpp
	g++ -g -Wall simulation.cpp -o simulation

clean:
	rm simulation
	rm -rf simulation.dSYM