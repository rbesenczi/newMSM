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
#ifndef NEWMESHREG_DISCRETEGROUPCOSTFUNCTION_H
#define NEWMESHREG_DISCRETEGROUPCOSTFUNCTION_H

#include "DiscreteCostFunction.h"

namespace newmeshreg {

class DiscreteGroupCostFunction : public NonLinearSRegDiscreteCostFunction {

public:
    //---INITIALISATION---//
    void initialize(int numNodes, int numLabels, int numPairs, int numTriplets) override {
        m_num_nodes = numNodes;
        m_num_labels = numLabels;
        m_num_pairs = numPairs;
        m_num_triplets = numTriplets;
    }

    void set_meshes(const std::vector<newresampler::Mesh>& source, const newresampler::Mesh& GRID, int num) override {
        _ORIG_MESHES = source;
        VERTICES_PER_SUBJ = GRID.nvertices();
        TRIPLETS_PER_SUBJ = GRID.ntriangles();
        _CONTROLMESHES.resize(num,GRID);
        subcorr = 0.1 * num;
    }

    void reset_CPgrid(const newresampler::Mesh& grid, int num) override { _CONTROLMESHES[num] = grid; }
    void set_patch_data(const std::vector<std::map<int,std::vector<double>>>& patches) override { patch_data = patches; }
    void set_masks(const newresampler::Mesh& m) override { _MASK = m; is_masked = true; }

    //---Compute costs---//
    double computeTripletCost(int triplet, int labelA, int labelB, int labelC) override;
    double computePairwiseCost(int pair, int labelA, int labelB) override;

private:
    std::vector<newresampler::Mesh> _CONTROLMESHES;
    std::vector<newresampler::Mesh> _ORIG_MESHES;
    newresampler::Mesh _MASK;

    std::vector<std::map<int,std::vector<double>>> patch_data;

    int TRIPLETS_PER_SUBJ = 0;
    int VERTICES_PER_SUBJ = 0;
    bool is_masked = false;
    double subcorr = 0.0;
};

} //namespace newmeshreg

#endif //NEWMESHREG_DISCRETEGROUPCOSTFUNCTION_H
