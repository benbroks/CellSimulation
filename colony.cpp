#include "colony.h"

Colony::Colony(int N, int X, int P, double SMin, double SMax, double R, double OR, double E, double M, double C, int binSize[], double binFlipRates[], bool V) {
    numCells = N;
    numBins = 51;
    numGenomes = 27634;
    minFlipRate = SMin;
    maxFlipRate = SMax;
    replaceRate = R;
    orderedReplaceRate = OR;
    expansionRate = E;
    verbose = V;
    neoplasticCycle = X;
    maxExpansionProportion = M;
    survivalRate = C;
    statFrequency = P;

    if (binFlipRates[0] == -1) {
        // Only perform exponential if first index is -1
        for(int i = 0; i < 51; i++) {
            if (i <= 25) {
                binFlipRates[i] = minFlipRate * pow(maxFlipRate / minFlipRate, double(i)/25);
            } else {
                binFlipRates[i] = binFlipRates[50-i];
            }
        }
    }

    for(int i = 0; i < numCells; i++) {
        healthyCells.insert(i);
    }
   
    replacePerTransition = OR * double(N);
    orderedReplacementCounter = 0;

    Cells = new Cell[N];
    if (verbose) {
        cout << "Begin Simulation." << endl;
    }
    for (int i = 0; i < N; i++) {
        Cells[i].setBinSize(binSize);
        Cells[i].generateGenome(binFlipRates);
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
    // Initialize
    for(int i = 0; i < numGenomes; i++) {
        avg[i] = 0;
    }
    // Find Overall Average
    for(int i = 0; i < numCells; i++) {
        for(int j = 0; j < numGenomes; j++) {
            avg[j] += Cells[i].getCpG(j);
        }
    }
    for (int i = 0; i < numGenomes; i++) {
        avg[i] = avg[i] / (2*numCells);
    }
}

// Calculate mean array for neoplastic cells
void Colony::findNeoplasticArray(double nAvg[]) {
    set <int> :: iterator itr1;
    // Initialize
    for(int i = 0; i < numGenomes; i++) {
        nAvg[i] = 0;
    }
    // Find Overall Average
    for (itr1 = neoplasticCells.begin(); itr1 != neoplasticCells.end(); itr1++) {
        for(int j = 0; j < numGenomes; j++) {
            nAvg[j] += Cells[*itr1].getCpG(j);
        }
    }
    // Normalize
    for (int i = 0; i < numGenomes; i++) {
        nAvg[i] = nAvg[i] / (2*neoplasticCells.size());
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

// Calculate Final Average Age
double Colony::findMeanAge() {
    double s = 0;
    for(int i = 0; i < numCells; i++) {
        s += Cells[i].getAge();
    }
    return s / numCells;
}

void Colony::transition(int T, string s_o_fp, string m_o_fp) {
    chrono::time_point<chrono::system_clock> start, end; 
    chrono::duration<double> elapsed_seconds;
    float r1;
    if (verbose) {
        start = chrono::system_clock::now(); 
    }
    for(int i = 0; i < T; i++) {
        // Print intermediate statistics to csv files
        if((statFrequency != -1) && (i != 0) && (i % statFrequency == 0)) {
            printStats(s_o_fp + "_" + to_string(i),i);
        }
        // Row/Cell Expansion
        if(i >= neoplasticCycle) {
            cellExpansion();
        } 
        // Ordered Cell Replacement
        if(replacePerTransition > 0) {
            int replaced = 0;
            while (replaced < replacePerTransition) {
                Cells[orderedReplacementCounter].cellReplacement();
                healthyCells.insert(orderedReplacementCounter);
                neoplasticCells.erase(orderedReplacementCounter);
                if (orderedReplacementCounter == (numCells-1)) {
                    orderedReplacementCounter = 0;
                } else {
                    orderedReplacementCounter++;
                }
                replaced++;
            }
        } 
        // Survival Replacement
        if (survivalRate < 1) {
            cellDeathAndReplacement();
        }

        // Normal Transition
        for(int j = 0; j < numCells; j++) {
            r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
            // Random Replacement
            if (r1 < replaceRate) {
                Cells[j].cellReplacement();
                healthyCells.insert(j);
                neoplasticCells.erase(j);
            }
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

void Colony::cellExpansion() {
    float r1;
    if(neoplasticCells.empty()) {
        int toAdd = rand() % numCells;
        neoplasticCells.insert(toAdd);
        healthyCells.erase(toAdd);
    } else {
        set <int> :: iterator itr1, itr2;
        for(itr1 = neoplasticCells.begin(); itr1 != neoplasticCells.end(); itr1++) {
            // Only expand if we are allowed to
            if (neoplasticCells.size() < maxExpansionProportion * numCells) {
                r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
                // Only expands if cell is randomly chosen
                if (r1 < expansionRate) {
                    // Picks random healthy cell to replace
                    int newCell = rand() % healthyCells.size();
                    itr2 = healthyCells.begin();
                    for(int i = 0; i < newCell; i++) {
                        itr2++;
                    }
                    Cells[*itr2] = Cells[*itr1];
                    Cells[*itr2].clearAge();
                    neoplasticCells.insert(*itr2);
                    healthyCells.erase(*itr2);
                }
            } else {
                return;
            }
        }
    }
}

// Very similar to cell expansion function
void Colony::cellDeathAndReplacement() {
    float r1;
    int newCell;
    set<int> liveCells, replacedCells;
    set <int> :: iterator itr1, itr2;
    // Choosing which cells survive vs. get replaced
    for (int i = 0; i < numCells; i++) {
        r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        // Survives
        if (r1 < survivalRate) {
            liveCells.insert(i);
        } else {
            // Replaced
            replacedCells.insert(i);
        }
    }

    // Randomly pick cell to replace old ones with
    for(itr1 = replacedCells.begin(); itr1 != replacedCells.end(); itr1++) {
        newCell = rand() % liveCells.size();
        itr2 = liveCells.begin();
        for(int i = 0; i < newCell; i++) {
            itr2++;
        }
        Cells[*itr2] = Cells[*itr1];
        Cells[*itr2].clearAge();
        neoplasticCells.insert(*itr2);
        healthyCells.erase(*itr2);
    }

}

void Colony::printStats(string o_fp, int numTransitions) {
    ofstream myfile;
    myfile.open(o_fp + ".csv"); 
    double avg[numGenomes];
    // Mean + Variance for All Cells
    findMeanArray(avg);
    double mu = findMean(avg);
    double var = findVariance(mu,avg);
    // Mean + Variance for Neoplastic Cells
    double nAvg[numGenomes];
    findNeoplasticArray(nAvg);
    double nMu = findMean(nAvg);
    double nVar = findVariance(nMu,nAvg);
    // Mean Age of Cells
    double ageMu = findMeanAge();
    // Printing to cout if need be
    if ((verbose) && (numTransitions == -1)) {
        cout << "Mean: " << mu << endl;
        cout << "Variance: " << var << endl;
        cout << "Mean Age: " << ageMu << endl;
        cout << "Neoplastic Cells: " << neoplasticCells.size() << endl;
        // We get nan if we don't have the zero check
        if (neoplasticCells.size() != 0) {
            cout << "Neoplastic Cell Mean: " << nMu << endl;
            cout << "Neoplastic Cell Variance: " << nVar << endl;
        }
    }
    if (numTransitions != -1) {
        myfile << "Completed Transitions," << numTransitions << endl;
    } else {
        myfile << "Final State" << endl;
    }
    myfile << "CgP Site,CgP Site Average for All Cells,CgP Site Average for Neoplastic Cells,,Mean," << mu << endl;
    for(int i = 0; i < numGenomes; i++) {
        myfile << i << "," << avg[i] << "," << nAvg[i];
        if (i == 0) {
            myfile << ",,Variance," << var;
        }
        else if (i == 1) {
            myfile << ",,Mean Age," << ageMu;
        } 
        else if (i == 2) {
            myfile << ",,Neoplastic Cells at Simulation End," << neoplasticCells.size();
        }
        else if (i == 3) {
            if (neoplasticCells.size() > 0) {
                myfile << ",,Neoplastic Cell Mean," << nMu;
            }
        }
        else if (i == 4) {
            if (neoplasticCells.size() > 0) {
                myfile << ",,Neoplastic Cell Variance," << nVar;
            }
        }
        myfile << endl;
    }
    myfile.close();
}

void Colony::printState(string o_fp, int numTransitions) {
    ofstream myfile;
    pair<int,int> p;
    // Print first half of CpG sites to one file
    myfile.open(o_fp + "1.csv");
    // Header
    if (numTransitions != -1) {
        myfile << "Completed Transitions," << numTransitions << endl;
    } else {
        myfile << "Final State" << endl;
    }
    myfile << "Cell \\ CpG Site,Age,";
    for (int i = 0; i < numGenomes/2; i++) {
        myfile << i << ",";
    }
    myfile << endl;
    // Print Cells
    for (int i = 0; i < numCells; i++) {
        myfile << i << "," << Cells[i].getAge() << ",";
        for(int j = 0; j < numGenomes / 2; j++) {
            myfile << Cells[i].getCpG(j) / float(2) << ",";
        }
        myfile << endl;
    }
    myfile.close();

    // Print back half of CpG sites to another file
    myfile.open(o_fp + "2.csv");
    // Header
    if (numTransitions != -1) {
        myfile << "Completed Transitions," << numTransitions << endl;
    } else {
        myfile << "Final State" << endl;
    }
    myfile << "Cell \\ CpG Site,Age,";
    for (int i = numGenomes/2; i < numGenomes; i++) {
        myfile << i << ",";
    }
    myfile << endl;
    // Print Cells
    for (int i = 0; i < numCells; i++) {
        myfile << i << "," << Cells[i].getAge() << ",";
        for(int j = numGenomes/2; j < numGenomes; j++) {

            myfile << Cells[i].getCpG(j) / float(2) << ",";
        }
        myfile << endl;
    }
    myfile.close();
}

