#include <fstream>
#include <utility>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <cmath>

#include "colony.h"


using namespace std;

string paramFile = "param.h";
string line, i_fp, m_o_fp, s_o_fp;
int N, T;
double S, R, OR, X, E, M;

// Goes through parameter file to obtain N,T,S,R,E values and relevant file paths

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
        else if (line[0] == 'O') {
            OR = stof(line.substr(3,line.length()));
        }
        else if (line[0] == 'X') {
            X = stof(line.substr(2,line.length()));
        }
        else if (line[0] == 'E') {
            E = stof(line.substr(2,line.length()));
        }
        else if (line[0] == 'M') {
            M = stof(line.substr(2,line.length()));
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
    ifstream input(i_fp);
    int i = 0;
    while(getline(input, line)) {
        if((line[0] >= '0') && (line[0] <= '9')) {
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
    }
    return binSize;
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

int main(int argc, char *argv[]){
    findParam();
    int binSize[51];
    findCPG(binSize);
    srand (static_cast <unsigned> (time(0)));
    // Instantiate
    Colony c = Colony(N, X, S, R, OR, E, M, binSize, true);
    // Transition
    c.transition(T);
    // Print Final Matrix and Statistics.
    c.printStats(s_o_fp);
    c.printFinalState(m_o_fp);
	return 1;
}