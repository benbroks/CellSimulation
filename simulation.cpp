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
int N, T, P, s;
double C, SA, SB, R, OR, X, E, M;

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
        else if (line[0] == 'C') {
            C = stof(line.substr(2,line.length()));
        }
        else if (line[0] == 'S') {
            if(line[1] == 'A') {
                SA = stof(line.substr(3,line.length()));
            } else {
                SB = stof(line.substr(3,line.length()));
            }
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
        else if (line[0] == 'A') {
            s = stoi(line.substr(2,line.length()));
        }
        else if (line[0] == 'P') {
            P = stoi(line.substr(2,line.length()));
            // Otherwise we'll get an ugly divide by zero error
            if (P == 0) {
                P = -1;
            }
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

int main(int argc, char *argv[]){
    findParam();
    int binSize[51];
    findCPG(binSize);
    if (s == -1) {
        srand (static_cast <unsigned> (time(0)));
    } else {
        srand(s);
    }
    // Instantiate
    Colony c = Colony(N, X, P, C, SA, SB, R, OR, E, M, binSize, true);
    // Transition
    c.transition(T, s_o_fp, m_o_fp);
    // Print Final Matrix and Statistics. -1 Indicates completed simulation.
    c.printStats(s_o_fp,-1);
    c.printState(m_o_fp,-1);
	return 1;
}