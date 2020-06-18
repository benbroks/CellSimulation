#include <vector>
#include <iostream>
#include <fstream>
#include <random>

using namespace std;

class Cell {
    private:
        int CpGBoxes, age;
        float CpGProportion,AflipRate,BflipRate;
        int bins[51];
        int findBin(int CpGSite);
        void randomCpGReplacement();
        char Genomes[27634];
    public:
        Cell();
        void generateGenome(float C, float SA, float SB);
        void cellReplacement();
        void transition();
        void setBinSize(int binSize[]);
        int getAge();
        char getCpG(int i);
        void clearAge();
        
};