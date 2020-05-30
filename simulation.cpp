#include <fstream>
#include <utility>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <cmath>
#include <chrono> 
#include "cell.h"


using namespace std;

string paramFile = "param.h";
string cpgFile = "cpg.csv";
string line, i_fp, m_o_fp, s_o_fp;
int N, T;
double S, R, E;

// Goes through parameter file to obtain N,T,S,R,E values.
// We can ignore the 'i' and 'o' values for now.

void findParam() {
    ifstream input(paramFile);
    while(getline(input, line)) {
        if (line[0] == 'N') {
            N = stoi(line.substr(2,line.length()));
        }
        else if (line[0] == 'T') {
            T = stoi(line.substr(2,line.length()));
        }
        else if (line[0] == 'S') {
            S = stof(line.substr(2,line.length()));
        }
        else if (line[0] == 'R') {
            R = stof(line.substr(2,line.length()));
        }
        else if (line[0] == 'E') {
            E = stof(line.substr(2,line.length()));
        }
        else if (line[0] == 'i') {
            int afterEqual = line.find("=") + 1;
            i_fp = line.substr(afterEqual,line.length());
        }
        else if (line[0] == 'm') {
            int afterEqual = line.find("=") + 1;
            m_o_fp = line.substr(afterEqual,line.length());
        }
        else if (line[0] == 's') {
            int afterEqual = line.find("=") + 1;
            s_o_fp = line.substr(afterEqual,line.length());
        }
    }
}

// Reads Excel file to obtain CpG distribution
int * findCPG(int * binSize) {
    ifstream input(cpgFile);
    int i = 0;
    while(getline(input, line)) {
        if(line[0] == '%') {
            continue;
        }
        int commaCount = 0;
        int j = 0;
        string sizeOfBin = "";
        while(commaCount < 5) {
            if(line[j] == ',') {
                commaCount ++;
            } else if (commaCount == 4) {
                sizeOfBin = sizeOfBin + line[j];
            }
            j ++;
        }
        binSize[i] = stoi(sizeOfBin);
        i ++;
    }
    return binSize;
}

// Calculate final mean
double findMean(Cell * Cells, ofstream & o) {
    Cell c;
    // Instantiate Average Array
    float avg[27634];
    for(int i = 0; i < 27634; i++) {
        avg[i] = 0;
    }
    // Find Overall Average
    for(int i = 0; i < N; i++) {
        c = Cells[i];
        for(int j = 0; j < 27634; j++) {
            avg[j] += Cells[i].getPair(j).first + Cells[i].getPair(j).second;
        }
    }
    float s = 0;
    o << "CgP Averages,";
    for(int i = 0; i < 27634;i++) {
        s += avg[i];
        o << avg[i] / (2*N) << ",";
    }
    o << endl;
    double mu = s / (double(2) * 27634 * N);
    o << "Mean," << mu << endl;
    cout << "Mean: " << mu << endl;
    return mu;
}

// Calculate Final Variance
void findVariance(Cell * Cells, double mean, ofstream & o) {
    double col_mean = 0;
    double var = 0;
    for(int i = 0; i < 27634; i++) {
        col_mean = 0;
        for(int j = 0; j < N; j++) {
            col_mean += (Cells[j].getPair(i).first + Cells[j].getPair(i).second) / double(2);
        }
        var += (col_mean/N - mean) * (col_mean/N - mean);
    }
    o << "Variance," << var / 27634 << endl;
    cout << "Variance: " << var / 27634 << endl;
}

// Calculates final CgP Concentration
void histogram(Cell * Cells) {
    int * Density = new int [N];
    int H[51];

    // Instantiating Occurence Arrays
    for (int i = 0; i < N; i++) {
        Density[i] = 0;
        if (i < 51) {
            H[i] = 0;
        }
    }
    // Calculating Sum of CgP Sites
    for (int i = 0; i < N; i++) {
        for(int j = 0; j < 27634; j++) {
            Density[i] += Cells[i].getPair(j).first + Cells[i].getPair(j).second;
        }
    }
    // Placing GcP Sites in Bins
    for (int i = 0; i < N; i++) {
        int param = int(double(Density[i]) / (2*27634) / 0.02);
        H[param] += 1;
    }
    // Printing Bins
    for (int i = 0; i < 51; i++) {
        cout << "Bin " << i << ": " << H[i] << endl;
    }
    delete [] Density;
}

void printFinalState(string o_fp, Cell * Cells) {
    ofstream myfile;
    myfile.open(o_fp); 
    // Print Header
    myfile << "Cell \\ CpG Site,";
    for (int i = 0; i < 27634; i++) {
        myfile << i << ",";
    }
    myfile << endl;
    // Print Cells
    for (int i = 0; i < N; i++) {
        myfile << i << ",";
        Cells[i].print(myfile);
    }
    myfile.close();
}

// Prints relevant statistics
void printStats(string o_fp, Cell * Cells) {
    ofstream myfile;
    myfile.open(o_fp); 
    // Print Header
    myfile << "CgP Sites,";
    for(int i = 0; i < 27634; i++) {
        myfile << i << ",";
    }
    myfile << endl;
    double mu = findMean(Cells,myfile);
    findVariance(Cells, mu,myfile);
}

int main(int argc, char *argv[]){
    findParam();
    int binSize[51];
    findCPG(binSize);
    srand (static_cast <unsigned> (time(0)));
    cout << "Begin Simulation." << endl;
    Cell * Cells = new Cell [N];
    for (int i = 0; i < N; i++) {
        Cells[i].setBinSize(binSize);
        Cells[i].generateGenome(S,R,E);
    }
    cout << "Cells instantiated." << endl;
    // Simulate T Transitions
    chrono::time_point<chrono::system_clock> start, end; 
    chrono::duration<double> elapsed_seconds;
    start = chrono::system_clock::now(); 
    for(int i = 0; i < T; i++) {
        for(int j = 0; j < N; j++) {
             Cells[j].transition();
        }
        end = chrono::system_clock::now(); 
        elapsed_seconds = end - start; 
        cout << "Completed " << i + 1 << " of " << T << " transitions. " << "Approximately " << T/float(i+1) * elapsed_seconds.count() - elapsed_seconds.count() <<" seconds remaining." << "\r";
        cout.flush();
    }
    cout << endl;
    // Print Final Matrix and Statistics.
    printStats(s_o_fp, Cells);
    printFinalState(m_o_fp, Cells);

    delete [] Cells;
	return 1;
}