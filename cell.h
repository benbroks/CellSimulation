#include <vector>
#include <iostream>
#include <fstream>
#include <random>

using namespace std;

class Cell {
    private:
        int CpGBoxes, age;
        float flipRate;
        int bins[51];
        int findBin(int CpGSite);
        void randomCpGReplacement();
        int Genomes[27634];
    public:
        Cell();
        void generateGenome(float S);
        void cellReplacement();
        void transition();
        void setBinSize(int binSize[]);
        int getAge();
        int getCpG(int i);
        void clearAge();
        
};