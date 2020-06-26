#include <vector>
#include <iostream>
#include <fstream>
#include <random>

using namespace std;

class Cell {
    private:
        int CpGBoxes, age;
        int bins[51];
        double * binFlipRates;
        int findBin(int CpGSite);
        void randomCpGReplacement();
        char Genomes[27634];
    public:
        Cell();
        void generateGenome(double * bFR);
        void cellReplacement();
        void transition();
        void setBinSize(int binSize[]);
        int getAge();
        char getCpG(int i);
        void clearAge();
        
};