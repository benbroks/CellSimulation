#include "cell.h"


Cell::Cell() {
    CpGBoxes = 27634;
    age = 0;
}

void Cell::generateGenome(double * f) {
    int currentBin = 0;
    flipRates = f;
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
    float r1, sampleProb, left, right;
    for(int i = 0; i < CpGBoxes; i++) {
        // Flip GcP Values with Bin Error Probabilities
        // flipRates determines the flip rate for this specific CpG site
        sampleProb = findBin(i) * flipRates[i] * 0.02;
        if (findBin(i) < 25) {
            left = sampleProb;
            right = flipRates[i];
        } else if (findBin(i) == 25) {
            // Center Bins have more variance
            left = flipRates[i];
            right = flipRates[i];
        } else {
            left = flipRates[i];
            right = flipRates[i] - sampleProb;
        }
        r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        if (Genomes[i] < 1) {
            // Methy Value
            if (r1 < left) {
                Genomes[i] += 1;
            }
        } else {
            // Demethy Value
            if (r1 < right) {
                Genomes[i] -= 1;
            }
        }
        r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        if (Genomes[i] <= 1) {
            // Methy Value
            if (r1 < left) {
                Genomes[i] += 1;
            }
        } else {
            // Demethy Value
            if (r1 < right) {
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

