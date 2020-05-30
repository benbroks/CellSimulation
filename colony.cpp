#include "colony.h"

// Calculate final mean
double Colony::findMean(ofstream & o) {
    Cell c;
    // Instantiate Average Array
    float avg[numGenomes];
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
    float s = 0;
    o << "CgP Averages,";
    for(int i = 0; i < numGenomes;i++) {
        s += avg[i];
        o << avg[i] / (2*numCells) << ",";
    }
    o << endl;
    double mu = s / (double(2) * numGenomes * numCells);
    o << "Mean," << mu << endl;
    if (verbose) {
        // Print to command line
        cout << "Mean: " << mu << endl;
    }
    return mu;
}

// Calculate Final Variance
void Colony::findVariance(double mean, ofstream & o) {
    double col_mean = 0;
    double var = 0;
    for(int i = 0; i < numGenomes; i++) {
        col_mean = 0;
        for(int j = 0; j < numCells; j++) {
            col_mean += (Cells[j].getPair(i).first + Cells[j].getPair(i).second) / double(2);
        }
        var += (col_mean/numCells - mean) * (col_mean/numCells - mean);
    }
    o << "Variance," << var / numGenomes << endl;
    if (verbose) {
        // Print to command line
        cout << "Variance: " << var / numGenomes << endl;
    }
}

Colony::~Colony() 
{ 
    delete [] Cells;
} 

Colony::Colony(int N, double S, double R, double OR, double E, int binSize[], bool V) {
    numCells = N;
    numBins = 51;
    numGenomes = 27634;
    flipRate = S;
    replaceRate = R;
    orderedReplaceRate = OR;
    expansionRate = E;
    verbose = V;
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

void Colony::transition(int T) {
    chrono::time_point<chrono::system_clock> start, end; 
    chrono::duration<double> elapsed_seconds;
    if (verbose) {
        start = chrono::system_clock::now(); 
    }
    for(int i = 0; i < T; i++) {
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
    // Print Header
    myfile << "CgP Sites,";
    for(int i = 0; i < numGenomes; i++) {
        myfile << i << ",";
    }
    myfile << endl;
    double mu = findMean(myfile);
    findVariance(mu,myfile);
    myfile.close();
}

void Colony::printFinalState(string o_fp) {
    ofstream myfile;
    myfile.open(o_fp); 
    // Print Header
    myfile << "Cell \\ CpG Site,";
    for (int i = 0; i < numGenomes; i++) {
        myfile << i << ",";
    }
    myfile << endl;
    // Print Cells
    for (int i = 0; i < numCells; i++) {
        myfile << i << ",";
        Cells[i].print(myfile);
    }
    myfile.close();
}

