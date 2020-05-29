#include "cell.h"
#include <iostream>
#include <random>

Cell::Cell() {
    CpGBoxes = 27634;
    age = 0;
}

void Cell::generateGenome(float S, float R, float E) {
    Genomes = vector< pair<int,int> > (CpGBoxes, pair<int,int>(0,0));
    int currentBin = 0;
    for (int i = 0; i < CpGBoxes; i++) {
        currentBin = findBin(i);
        // Setting first column value
        float r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        if (r1 < currentBin * 0.02) {
            Genomes[i].first = 1;
        } else {
            Genomes[i].first = 0;
        }
        // Setting second column value
        r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        if (r1 < currentBin * 0.02) {
            Genomes[i].second = 1;
        } else {
            Genomes[i].second = 0;
        }
    }
}

void Cell::transition() {
    float r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    int currentBin = 0;
    if (r1 < R) {
        // Generate completely new Genome with Prob. R (Random Replacement)
        for (int i = 0; i < CpGBoxes; i++) {
            currentBin = findBin(i);
            // Setting first column value
            float r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
            if (r1 < currentBin * 0.02) {
                Genomes[i].first = 1;
            } else {
                Genomes[i].first = 0;
            }
            // Setting second column value
            r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
            if (r1 < currentBin * 0.02) {
                Genomes[i].second = 1;
            } else {
                Genomes[i].second = 0;
            }
        }
    } else {
        for(int i = 0; i < CpGBoxes; i++) {
            currentBin = findBin(i);
            // Flip GcP Values with Bin Error Probabilities
            r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
            if (Genomes[i].first == 0) {
                // Methy Value
                if (r1 < currentBin*S*0.02) {
                    Genomes[i].first = 1;
                }
            } else {
                // Demethy Value
                if (r1 < (1-currentBin*0.02)*S) {
                    Genomes[i].first = 0;
                }
            }
            r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
            if (Genomes[i].second == 0) {
                // Methy Value
                if (r1 < currentBin*S*0.02) {
                    Genomes[i].second = 1;
                }
            } else {
                // Demethy Value
                if (r1 < (1-currentBin*0.02)*S) {
                    Genomes[i].second = 0;
                }
            }
        }
    }

    
}

void Cell::setBinSize(int binSize[]) {
    for(int i = 0; i < 51; i++) {
        bins[i] = binSize[i];
    }
}

pair<int,int> Cell::getPair(int i) {
    return Genomes[i];
}

int Cell::findBin(int CpGSite) {
    int b = -1;
    while (CpGSite > 0) {
        b++;
        CpGSite -= bins[b];
    }
    return b;
}

