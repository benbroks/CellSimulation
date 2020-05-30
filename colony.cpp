#include "colony.h"

Colony::Colony(int N, double S, double R, double OR, double E, int binSize[], bool V) {
    numCells = N;
    numBins = 51;
    numGenomes = 27634;
    flipRate = S;
    replaceRate = R;
    orderedReplaceRate = OR;
    expansionRate = E;
    verbose = V;

    replacePerTransition = int(OR * double(N));
    orderedReplacementCounter = 0;
    if (replacePerTransition != 0) {
        if ((N % replacePerTransition) != 0) {
            cout << "Ordered Replacement Proportion does not evenly divide N. Undefined behavior to follow." << endl;
        }
    }

    Cells = new Cell[N];
    if (verbose) {
        cout << "Begin Simulation." << endl;
    }
    for (int i = 0; i < N; i++) {
        Cells[i].setBinSize(binSize);
        Cells[i].generateGenome(S,R,E);
    }
    if (verbose) {
        cout << "Cells instantiated." << endl;
    }
}

Colony::~Colony() 
{ 
    delete [] Cells;
} 

// Calculate mean array
void Colony::findMeanArray(double * avg) {
    Cell c;
    for(int i = 0; i < numGenomes; i++) {
        avg[i] = 0;
    }
    // Find Overall Average
    for(int i = 0; i < numCells; i++) {
        c = Cells[i];
        for(int j = 0; j < numGenomes; j++) {
            avg[j] += Cells[i].getPair(j).first + Cells[i].getPair(j).second;
        }
    }
    for (int i = 0; i < numGenomes; i++) {
        avg[i] = avg[i] / (2*numCells);
    }
}

// Calculate final mean
double Colony::findMean(double avg[]) {
    double s = 0;
    for(int i = 0; i < numGenomes; i++) {
        s += avg[i];
    }
    return s / (numGenomes);
}

// Calculate Final Variance
double Colony::findVariance(double mean, double avg[]) {
    double var = 0;
    for(int i = 0; i < numGenomes; i++) {
        var += (avg[i] - mean) * (avg[i] - mean);
    }
    return var / numGenomes;
}

void Colony::transition(int T) {
    chrono::time_point<chrono::system_clock> start, end; 
    chrono::duration<double> elapsed_seconds;
    if (verbose) {
        start = chrono::system_clock::now(); 
    }
    for(int i = 0; i < T; i++) {
        // Ordered Cell Replacement
        if(replacePerTransition > 0) {
            for(int j = orderedReplacementCounter; j < orderedReplacementCounter + replacePerTransition; j++) {
                Cells[j].cellReplacement();
            }
            orderedReplacementCounter = (orderedReplacementCounter + replacePerTransition) % numCells;
        } 
        // Normal Transition
        for(int j = 0; j < numCells; j++) {
             Cells[j].transition();
        }
        if (verbose) {
            end = chrono::system_clock::now(); 
            elapsed_seconds = end - start; 
            cout << "Completed " << i + 1 << " of " << T << " transitions. " << "Approximately " << T/float(i+1) * elapsed_seconds.count() - elapsed_seconds.count() <<" seconds remaining." << "\r";
            cout.flush();
        }
    }
    if(verbose) {
        cout << endl;
    }
}

void Colony::printStats(string o_fp) {
    ofstream myfile;
    myfile.open(o_fp); 
    double avg[numGenomes];
    findMeanArray(avg);
    double mu = findMean(avg);
    double var = findVariance(mu,avg);
    if (verbose) {
        cout << "Mean: " << mu << endl;
        cout << "Variance: " << var << endl;
    }
    myfile << "CgP Site,CgP Site Average,,Mean," << mu << endl;
    for(int i = 0; i < numGenomes; i++) {
        myfile << i << "," << avg[i];
        if (i == 0) {
            myfile << ",,Variance," << var;
        }
        myfile << endl;
    }
    myfile.close();
}

void Colony::printFinalState(string o_fp) {
    ofstream myfile;
    pair<int,int> p;
    // Print first half of CpG sites to one file
    myfile.open(o_fp + "1.csv");
    // Header
    myfile << "Cell \\ CpG Site,";
    for (int i = 0; i < numGenomes/2; i++) {
        myfile << i << ",";
    }
    myfile << endl;
    // Print Cells

    for (int i = 0; i < numCells; i++) {
        myfile << i << ",";
        for(int j = 0; j < numGenomes / 2; j++) {
            p = Cells[i].getPair(j);
            myfile << (p.first + p.second) / float(2) << ",";
        }
        myfile << endl;
    }
    myfile.close();

    // Print back half of CpG sites to another file
    myfile.open(o_fp + "2.csv");
    // Header
    myfile << "Cell \\ CpG Site,";
    for (int i = numGenomes/2; i < numGenomes; i++) {
        myfile << i << ",";
    }
    myfile << endl;
    // Print Cells
    for (int i = 0; i < numCells; i++) {
        myfile << i << ",";
        for(int j = numGenomes/2; j < numGenomes; j++) {
            p = Cells[i].getPair(j);
            myfile << (p.first + p.second) / float(2) << ",";
        }
        myfile << endl;
    }
    myfile.close();
}

