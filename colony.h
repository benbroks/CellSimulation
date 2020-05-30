#include <chrono> 
#include "cell.h"

using namespace std;

class Colony {
    private:
        int numCells, numGenomes, numBins, replacePerTransition, orderedReplacementCounter;
        double flipRate,replaceRate,orderedReplaceRate,expansionRate;
        bool verbose;
        Cell * Cells;
        void findMeanArray(double avg[]);
        double findMean(double avg[]);
        double findVariance(double mean, double avg[]);
    public:
        Colony(int N, double S, double R, double OR, double E, int binSize[], bool v);
        ~Colony();
        void transition(int T);
        void printStats(string o_fp);
        void printFinalState(string o_fp);
};