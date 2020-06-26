#include <chrono> 
#include <set>
#include <unordered_set>
#include <filesystem>
#include "cell.h"

using namespace std;

class Colony {
    private:
        int numCells, numGenomes, numBins, orderedReplacementCounter, neoplasticCycle, statFrequency;
        double minFlipRate, maxFlipRate, replaceRate,orderedReplaceRate,expansionRate,maxExpansionProportion,replacePerTransition;
        double * binFlipRates;
        bool verbose;
        set<int> neoplasticCells, healthyCells;
        Cell * Cells;
        void findMeanArray(double avg[]);
        double findMean(double avg[]);
        double findVariance(double mean, double avg[]);
        double findMeanAge();
        void findNeoplasticArray(double nAvg[]);
    public:
        Colony(int N, int X, int P, double SMin, double SMax, double R, double OR, double E, double M, int binSize[], bool v);
        ~Colony();
        void transition(int T, string s_o_fp, string m_o_fp);
        void cellExpansion();
        void printStats(string o_fp, int numTransitions);
        void printState(string o_fp, int numTransitions);
};