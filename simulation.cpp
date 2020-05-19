#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <cmath>
#include <set>

using namespace std;

string paramFile = "param.h";
string line, i_fp, o_fp;
int N, T;
float S, R, E;

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

int main(int argc, char *argv[]){
    findParam();
	return 1;
}