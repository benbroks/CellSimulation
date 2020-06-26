#include "cell.h"


Cell::Cell() {
    CpGBoxes = 27634;
    age = 0;
}

void Cell::generateGenome(double * bFR) {
    int currentBin = 0;
    binFlipRates = bFR;
    float r1;
    for (int i = 0; i < CpGBoxes; i++) {
        currentBin = findBin(i);
        r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        // p = 0.02 * currentBin
        if (r1 < 0.0004 * currentBin * currentBin) {
            // Prob of p*p
            Genomes[i] = 2;
        } else if (r1 < 0.04 * currentBin * (1 - .01 * currentBin)) {
            // Prob of 2*(1-p)*p
            Genomes[i] = 1;
        } else {
            // Prob of (1-p)*(1-p)
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
        // p = 0.02 * currentBin
        if (r1 < 0.0004 * currentBin * currentBin) {
            // Prob of p*p
            Genomes[i] = 2;
        } else if (r1 < 0.04 * currentBin * (1 - .01 * currentBin)) {
            // Prob of 2*(1-p)*p
            Genomes[i] = 1;
        } else {
            // Prob of (1-p)*(1-p)
            Genomes[i] = 0;
        }
    }
    age = 0;
}

void Cell::randomCpGReplacement() {
    float r1;
    double methy,demethy;
    int bin;
    for(int i = 0; i < CpGBoxes; i++) {
        // Flip GcP Values with Bin Error Probabilities
        
        // S is higher for "middle" bin sizes
        bin = findBin(i);

        methy = binFlipRates[i] * double(bin) / 50;
        demethy = binFlipRates[i] * (1 - double(bin) / 50);
        
        // One side of CpG Site
        r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        if (Genomes[i] < 1) {
            // Methy Value
            if (r1 < methy) {
                Genomes[i] += 1;
            }
        } else {
            // Demethy Value
            if (r1 < demethy) {
                Genomes[i] -= 1;
            }
        }
        // Other side of CpG Site
        r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        if (Genomes[i] <= 1) {
            // Methy Value
            if (r1 < methy) {
                Genomes[i] += 1;
            }
        } else {
            // Demethy Value
            if (r1 < demethy) {
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

char Cell::getCpG(int i) {
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

