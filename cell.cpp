#include "cell.h"


Cell::Cell() {
    CpGBoxes = 27634;
    age = 0;
}

void Cell::generateGenome(float C, float SA, float SB) {
    int currentBin = 0;
    CpGProportion = C;
    AflipRate = SA;
    BflipRate = SB;
    float r1;
    for (int i = 0; i < CpGBoxes; i++) {
        currentBin = abs(findBin(i)) - 1;
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
        currentBin = abs(findBin(i)) - 1;
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
    float r1, sampleProb;
    // Indicates which group we're in [A or B]
    bool cpgSide;
    for(int i = 0; i < CpGBoxes; i++) {
        // Flip GcP Values with Bin Error Probabilities
        // Checks to see which group the site falls into - this determines the error rate
        sampleProb = findBin(i);
        if (sampleProb < 0) {
            sampleProb += 1;
            sampleProb *= -1 * AflipRate;
            cpgSide = true;
        } else {
            sampleProb -= 1;
            sampleProb *= BflipRate;
            cpgSide = false;
        }
        sampleProb *= 0.02;

        r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        if (Genomes[i] < 1) {
            // Methy Value
            if (r1 < sampleProb) {
                Genomes[i] += 1;
            }
        } else {
            // Demethy Value
            if (cpgSide) {
                if (r1 < AflipRate - sampleProb) {
                    Genomes[i] -= 1;
                }
            } else {
                if (r1 < BflipRate - sampleProb) {
                    Genomes[i] -= 1;
                }
            }
        }
        r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        if (Genomes[i] <= 1) {
            // Methy Value
            if (r1 < sampleProb) {
                Genomes[i] += 1;
            }
        } else {
            // Demethy Value
            if (cpgSide) {
                if (r1 < AflipRate - sampleProb) {
                    Genomes[i] -= 1;
                }
            } else {
                if (r1 < BflipRate - sampleProb) {
                    Genomes[i] -= 1;
                }
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
    // Added a slight "hack" here 
    // First, we check how many CpG sites deep the cell is into a specific bin
    // This is used to determine whether the site falls into the A or B group
    // We return a negative number if in group A and positive if in group B
    // However, 0 can't be turned into a negative number so we have to add 1 first
    // The 1 is subtracted after returning the absolute value of the result if the group isn't needed
    if ((-1 * CpGSite) <= int(CpGProportion * bins[b])) {
        return -1 * (b+1);
    }
    return b+1;
}

void Cell::clearAge() {
    age = 0;
}

