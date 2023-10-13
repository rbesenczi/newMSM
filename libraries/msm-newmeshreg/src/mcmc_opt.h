/*
Copyright (c) 2023 King's College London, MeTrICS Lab, Renato Besenczi

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#ifndef NEWMESHREG_MCMC_OPT_H
#define NEWMESHREG_MCMC_OPT_H

#include "DiscreteModel.h"

namespace newmeshreg {

class MCMC {
public:
    static double optimise(const std::shared_ptr<NonLinearSRegDiscreteModel>& energy, bool verbose, int mciters) {

        int* labeling = energy->getLabeling();
        const int num_nodes = energy->getNumNodes();
        const double* unary_costs = energy->getCostFunction()->getUnaryCosts();
        const int* triplets = energy->getTriplets();

        if(verbose) { std::cout << "Initial "; energy->evaluateTotalCostSum(); }

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distribution(0, energy->getNumLabels()-1);

        for(int i = 0; i < mciters; ++i)
        {
            for(int triplet = 0; triplet < energy->getNumTriplets(); ++triplet)
            {
                int label = distribution(gen);

                const int nodeA = triplets[triplet*3];
                const int nodeB = triplets[triplet*3+1];
                const int nodeC = triplets[triplet*3+2];

                if ((energy->computeTripletCost(triplet, label, labeling[nodeB], labeling[nodeC]) +
                    unary_costs[label * num_nodes + nodeA]) <
                    (energy->computeTripletCost(triplet, labeling[nodeA], labeling[nodeB], labeling[nodeC]) +
                    unary_costs[labeling[nodeA] * num_nodes + nodeA]))
                        labeling[nodeA] = label;
            }

            if (verbose && i % 1000 == 0 && i > 0)
            {
                std::cout << "MC iter " << i << '/' << mciters << '\t';
                energy->evaluateTotalCostSum();
            }
        }

        if (verbose) std::cout << "MC iter " << mciters << '/' << mciters << '\t';

        return energy->evaluateTotalCostSum();
    }
};

} //namespace newmeshreg

#endif //NEWMESHREG_MCMC_OPT_H
