#include <chrono> 
#include "cell.h"

using namespace std;

class Colony {
    private:
        int numCells, numGenomes, numBins;
        double flipRate,replaceRate,orderedReplaceRate,expansionRate;
        bool verbose;
        Cell * Cells;
        double findMean(ofstream & o);
        void findVariance(double mean, ofstream & o);
    public:
        Colony(int N, double S, double R, double OR, double E, int binSize[], bool v);
        ~Colony();
        void transition(int T);
        void printStats(string o_fp);
        void printFinalState(string o_fp);
};