#include <vector>

using namespace std;

class Cell {
    private:
        int CpGBoxes, age;
        float S, R, E;
        int bins[51];
        vector< pair<int, int> > Genomes;
    public:
        Cell();
        void generateGenome(float S, float R, float E);
        void transition();
        void setBinSize(int binSize[]);
        pair<int,int> getPair(int i);
        
};