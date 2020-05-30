#include "cell.h"


Cell::Cell() {
    CpGBoxes = 27634;
    age = 0;
}

void Cell::generateGenome(float S) {
    int currentBin = 0;
    flipRate = S;
    float r1;
    for (int i = 0; i < CpGBoxes; i++) {
        currentBin = findBin(i);
        r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        if (r1 < 0.0004 * currentBin * currentBin) {
            Genomes[i] = 2;
        
        } else if (r1 < 0.04 * currentBin * (1 - .01 * currentBin)) {
            Genomes[i] = 1;
        } else {
            Genomes[i] = 0;
        }
    }
    age = 0;
}

void Cell::cellReplacement() {
    // Generate completely new Genome
    int currentBin = 0;
    float r1;
    for (int i = 0; i < CpGBoxes; i++) {
        currentBin = findBin(i);
        r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        if (r1 < 0.0004 * currentBin * currentBin) {
            Genomes[i] = 2;
        
        } else if (r1 < 0.04 * currentBin * (1 - .01 * currentBin)) {
            Genomes[i] = 1;
        } else {
            Genomes[i] = 0;
        }
    }
    age = 0;
}

void Cell::randomCpGReplacement() {
    int currentBin = 0;
    float r1;
    for(int i = 0; i < CpGBoxes; i++) {
        currentBin = findBin(i);
        // Flip GcP Values with Bin Error Probabilities
        r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        if (Genomes[i] < 1) {
            // Methy Value
            if (r1 < currentBin*flipRate*0.02) {
                Genomes[i] += 1;
            }
        } else {
            // Demethy Value
            if (r1 < (1-currentBin*0.02)*flipRate) {
                Genomes[i] -= 1;
            }
        }
        r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        if (Genomes[i] <= 1) {
            // Methy Value
            if (r1 < currentBin*flipRate*0.02) {
                Genomes[i] += 1;
            }
        } else {
            // Demethy Value
            if (r1 < (1-currentBin*0.02)*flipRate) {
                Genomes[i] -= 1;
            }
        }
    }
}

void Cell::transition() {
    randomCpGReplacement();
    age++;
}

void Cell::setBinSize(int binSize[]) {
    for(int i = 0; i < 51; i++) {
        bins[i] = binSize[i];
    }
}

int Cell::getAge() {
    return age;
}

short Cell::getCpG(int i) {
    return Genomes[i];
}

int Cell::findBin(int CpGSite) {
    int b = -1;
    while (CpGSite >= 0) {
        b++;
        CpGSite -= bins[b];
    }
    return b;
}

void Cell::clearAge() {
    age = 0;
}

