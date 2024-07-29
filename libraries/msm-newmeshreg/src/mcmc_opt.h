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
        const int num_triplets = energy->getNumTriplets();
        const int num_labels = energy->getNumLabels();
        const double* unary_costs = energy->getCostFunction()->getUnaryCosts();
        const double* triplet_costs = energy->getCostFunction()->getTripletCosts(); //TODO pre-calculate triplet costs
        const int* triplets = energy->getTriplets();

        if(verbose) { std::cout << "Initial "; energy->evaluateTotalCostSum(); }

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distribution(0, energy->getNumLabels()-1);

        for(int i = 0; i < mciters; ++i)
        {
            for(int triplet = 0; triplet < energy->getNumTriplets(); ++triplet) {
                int label = distribution(gen);
                std::vector<double> costs(8,0.0);

                const int nodeA = triplets[triplet*3];
                const int nodeB = triplets[triplet*3+1];
                const int nodeC = triplets[triplet*3+2];
                // TODO try all combinations -- move all vertices of the triplet to the label (see Fusion.h)
                /*
                if ((energy->computeTripletCost(triplet, label, labeling[nodeB], labeling[nodeC]) +
                    unary_costs[label * num_nodes + nodeA]) <
                    (energy->computeTripletCost(triplet, labeling[nodeA], labeling[nodeB], labeling[nodeC]) +
                    unary_costs[labeling[nodeA] * num_nodes + nodeA]))
                        labeling[nodeA] = label;
                */
                costs[0] = triplet_costs[triplet * num_triplets + labeling[nodeA] * num_labels + labeling[nodeB] * num_labels + labeling[nodeC]]
                           + unary_costs[labeling[nodeA] * num_nodes + nodeA]
                           + unary_costs[labeling[nodeB] * num_nodes + nodeB]
                           + unary_costs[labeling[nodeC] * num_nodes + nodeC];
                costs[1] = triplet_costs[triplet * num_triplets + labeling[nodeA] * num_labels + labeling[nodeB] * num_labels + label]
                           + unary_costs[labeling[nodeA] * num_nodes + nodeA]
                           + unary_costs[labeling[nodeB] * num_nodes + nodeB]
                           + unary_costs[label * num_nodes + nodeC];
                costs[2] = triplet_costs[triplet * num_triplets + labeling[nodeA] * num_labels + label * num_labels + labeling[nodeC]]
                           + unary_costs[labeling[nodeA] * num_nodes + nodeA]
                           + unary_costs[label * num_nodes + nodeB]
                           + unary_costs[labeling[nodeC] * num_nodes + nodeC];
                costs[3] = triplet_costs[triplet * num_triplets + labeling[nodeA] * num_labels + label * num_labels + label]
                            + unary_costs[labeling[nodeA] * num_nodes + nodeA]
                            + unary_costs[label * num_nodes + nodeB]
                            + unary_costs[label * num_nodes + nodeC];
                costs[4] = triplet_costs[triplet * num_triplets + label * num_labels + labeling[nodeB] * num_labels + labeling[nodeC]]
                           + unary_costs[label * num_nodes + nodeA]
                           + unary_costs[labeling[nodeB] * num_nodes + nodeB]
                           + unary_costs[labeling[nodeC] * num_nodes + nodeC];
                costs[5] = triplet_costs[triplet * num_triplets + label * num_labels + labeling[nodeB] * num_labels + label]
                           + unary_costs[label * num_nodes + nodeA]
                           + unary_costs[labeling[nodeB] * num_nodes + nodeB]
                           + unary_costs[label * num_nodes + nodeC];
                costs[6] = triplet_costs[triplet * num_triplets + label * num_labels + label * num_labels + labeling[nodeC]]
                           + unary_costs[label * num_nodes + nodeA]
                           + unary_costs[label * num_nodes + nodeB]
                           + unary_costs[labeling[nodeC] * num_nodes + nodeC];
                costs[7] = triplet_costs[triplet * num_triplets + label * num_labels + label * num_labels + label]
                           + unary_costs[label * num_nodes + nodeA]
                           + unary_costs[label * num_nodes + nodeB]
                           + unary_costs[label * num_nodes + nodeC];

                auto it = std::min_element(costs.begin(), costs.end());
                auto index = std::distance(costs.begin(), it);

                switch(index) {
                    case 0:
                        continue;
                    case 1:
                        labeling[nodeC] = label;
                        break;
                    case 2:
                        labeling[nodeB] = label;
                        break;
                    case 3:
                        labeling[nodeB] = labeling[nodeC] = label;
                        break;
                    case 4:
                        labeling[nodeA] = label;
                        break;
                    case 5:
                        labeling[nodeA] = labeling[nodeC] = label;
                        break;
                    case 6:
                        labeling[nodeA] = labeling[nodeB] = label;
                        break;
                    case 7:
                        labeling[nodeA] = labeling[nodeB] = labeling[nodeC] = label;
                        break;
                    default:
                        std::cout << index << std::endl;
                        throw MeshregException("Unknown error");
                }

            }

            if (verbose && i % 10000 == 0 && i > 0)
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
