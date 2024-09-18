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
        const int num_labels = energy->getNumLabels();
        const int num_triplets = energy->getNumTriplets();
        const double dist_param = energy->getMCParam();
        const int* triplets = energy->getTriplets();
        const double* unary_costs = energy->getCostFunction()->getUnaryCosts();
        const auto& tcosts = energy->getCostFunction()->getTCosts();
        int label = 0;
        std::vector<double> costs(8,0.0);

        if(verbose) { std::cout << "Initial "; energy->evaluateTotalCostSum(); std::cout << "Running Monte Carlo simulation..." << std::endl; }

        std::random_device rd;
        std::mt19937 gen(rd());
        std::geometric_distribution<> distribution(dist_param);

        for(int i = 0; i < mciters; ++i) {
            for(int triplet = 0; triplet < num_triplets; ++triplet) {
                do { label = distribution(gen); } while(label >= num_labels);

                const int nodeA = triplets[triplet*3];
                const int nodeB = triplets[triplet*3+1];
                const int nodeC = triplets[triplet*3+2];

                costs[0] = tcosts[triplet][labeling[nodeA]][labeling[nodeB]][labeling[nodeC]]
                           + (unary_costs[labeling[nodeA] * num_nodes + nodeA]
                           + unary_costs[labeling[nodeB] * num_nodes + nodeB]
                           + unary_costs[labeling[nodeC] * num_nodes + nodeC])/3.0;
                costs[1] = tcosts[triplet][labeling[nodeA]][labeling[nodeB]][label]
                           + (unary_costs[labeling[nodeA] * num_nodes + nodeA]
                           + unary_costs[labeling[nodeB] * num_nodes + nodeB]
                           + unary_costs[label * num_nodes + nodeC])/3.0;
                costs[2] = tcosts[triplet][labeling[nodeA]][label][labeling[nodeC]]
                           + (unary_costs[labeling[nodeA] * num_nodes + nodeA]
                           + unary_costs[label * num_nodes + nodeB]
                           + unary_costs[labeling[nodeC] * num_nodes + nodeC])/3.0;
                costs[3] = tcosts[triplet][labeling[nodeA]][label][label]
                           + (unary_costs[labeling[nodeA] * num_nodes + nodeA]
                           + unary_costs[label * num_nodes + nodeB]
                           + unary_costs[label * num_nodes + nodeC])/3.0;
                costs[4] = tcosts[triplet][label][labeling[nodeB]][labeling[nodeC]]
                           + (unary_costs[label * num_nodes + nodeA]
                           + unary_costs[labeling[nodeB] * num_nodes + nodeB]
                           + unary_costs[labeling[nodeC] * num_nodes + nodeC])/3.0;
                costs[5] = tcosts[triplet][label][labeling[nodeB]][label]
                           + (unary_costs[label * num_nodes + nodeA]
                           + unary_costs[labeling[nodeB] * num_nodes + nodeB]
                           + unary_costs[label * num_nodes + nodeC])/3.0;
                costs[6] = tcosts[triplet][label][label][labeling[nodeC]]
                           + (unary_costs[label * num_nodes + nodeA]
                           + unary_costs[label * num_nodes + nodeB]
                           + unary_costs[labeling[nodeC] * num_nodes + nodeC])/3.0;
                costs[7] = tcosts[triplet][label][label][label]
                           + (unary_costs[label * num_nodes + nodeA]
                           + unary_costs[label * num_nodes + nodeB]
                           + unary_costs[label * num_nodes + nodeC])/3.0;

                switch(std::distance(costs.begin(), std::min_element(costs.begin(), costs.end()))) {
                    case 0:
                        break;
                    case 1:
                        labeling[nodeC] = label;
                        break;
                    case 2:
                        labeling[nodeB] = label;
                        break;
                    case 3:
                        labeling[nodeB] = label;
                        labeling[nodeC] = label;
                        break;
                    case 4:
                        labeling[nodeA] = label;
                        break;
                    case 5:
                        labeling[nodeA] = label;
                        labeling[nodeC] = label;
                        break;
                    case 6:
                        labeling[nodeA] = label;
                        labeling[nodeB] = label;
                        break;
                    case 7:
                        labeling[nodeA] = label;
                        labeling[nodeB] = label;
                        labeling[nodeC] = label;
                        break;
                    default:
                        throw MeshregException("Unknown error");
                }
            }
            if (verbose && i % 10000 == 0) {
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
