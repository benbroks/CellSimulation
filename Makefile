all: simulation

simulation: simulation.cpp cell.cpp colony.cpp
	g++ -g -Wall simulation.cpp cell.cpp colony.cpp -o simulation

clean:
	rm simulation
	rm -rf simulation.dSYM
	rm matrix1.csv
	rm matrix2.csv
	rm stats.csv
	rm .DS_Store