#include <chrono> 
#include <set>
#include "cell.h"

using namespace std;

class Colony {
    private:
        int numCells, numGenomes, numBins, replacePerTransition, orderedReplacementCounter, neoplasticCycle;
        double flipRate,replaceRate,orderedReplaceRate,expansionRate,maxExpansionProportion;
        bool verbose;
        set<int> neoplasticCells, healthyCells;
        Cell * Cells;
        void findMeanArray(double avg[]);
        double findMean(double avg[]);
        double findVariance(double mean, double avg[]);
    public:
        Colony(int N, int X, double S, double R, double OR, double E, double M, int binSize[], bool v);
        ~Colony();
        void transition(int T);
        void cellExpansion();
        void printStats(string o_fp);
        void printFinalState(string o_fp);
};