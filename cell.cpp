#include "cell.h"
#include <iostream>
#include <random>

Cell::Cell() {
    CpGBoxes = 27634;
    age = 0;
}

void Cell::generateGenome(float S, float R, float E) {
    Genomes = vector< pair<int,int> > (CpGBoxes, pair<int,int>(0,0));
    for (int i = 0; i < CpGBoxes; i++) {
        float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        if (r > 0.75) {
            Genomes[i].first = 1;
            Genomes[i].second = 1;
        } else if (r > 0.5) {
            Genomes[i].first = 1;
            Genomes[i].second = 0;
        } else if (r > 0.25) {
            Genomes[i].first = 0;
            Genomes[i].second = 1;
        }
    }
}

pair<int,int> Cell::getPair(int i) {
    return Genomes[i];
}

void Cell::transition() {
    float r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    if (r1 < R) {
        // Generate completely new Genome with Prob. R (Random Replacement)
        for (int i = 0; i < CpGBoxes; i++) {
            float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
            if (r > 0.75) {
                Genomes[i].first = 1;
                Genomes[i].second = 1;
            } else if (r > 0.5) {
                Genomes[i].first = 1;
                Genomes[i].second = 0;
            } else if (r > 0.25) {
                Genomes[i].first = 0;
                Genomes[i].second = 1;
            } else {
                Genomes[i].first = 0;
                Genomes[i].second = 0;
            }
        }
    } else {
        for(int i = 0; i < CpGBoxes; i++) {
            // Flip GcP Values with Prob. S
            r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
            if(r1 < S) {
                Genomes[i].first = Genomes[i].first*(-1) + 1;
            }
            r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
            if(r1 < S) {
                Genomes[i].second = Genomes[i].second*(-1) + 1;
            }
        }
    }

    
}