#include <iostream>
#include <fstream>
#include <utility>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <cmath>
#include "cell.h"

using namespace std;

string paramFile = "param.h";
string line, i_fp, o_fp;
int N, T;
float S, R, E;
vector<float> binProb;

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
        else if (line[0] == 'o') {
            int afterEqual = line.find("=") + 1;
            o_fp = line.substr(afterEqual,line.length());
        }
    }
}

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
        int param = int(float(Density[i]) / (2*27634) / 0.02);
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
    srand (static_cast <unsigned> (time(0)));
    Cell * Cells = new Cell [N];
    for (int i = 0; i < N; i++) {
        Cells[i].generateGenome(S,R,E);
    }
    // Simulate T Transitions
    for(int i = 0; i < T; i++) {
        for(int j = 0; j < N; j++) {
             Cells[j].transition();
        }
    }
    histogram(Cells);

    delete [] Cells;
	return 1;
}