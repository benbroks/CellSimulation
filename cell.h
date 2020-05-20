#include <vector>

using namespace std;

class Cell {
    private:
        int CpGBoxes, age;
        float S, R, E;
        vector< pair<int, int> > Genomes;
    public:
        Cell();
        void generateGenome(float S, float R, float E);
        pair<int,int> getPair(int i);
        void transition();
};