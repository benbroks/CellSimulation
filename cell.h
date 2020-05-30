#include <vector>
#include <iostream>
#include <fstream>
#include <random>

using namespace std;

class Cell {
    private:
        int CpGBoxes, age;
        float S, R;
        int bins[51];
        int findBin(int CpGSite);
        void randomCpGReplacement();
        vector< pair<int, int> > Genomes;
    public:
        Cell();
        void generateGenome(float S, float R);
        void cellReplacement();
        void transition();
        void setBinSize(int binSize[]);
        pair<int,int> getPair(int i);
        int getAge();
        void clearAge();
        
};