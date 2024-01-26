#ifndef NEWMESHREG_DISCRETEGROUPCOFUNCTION_H
#define NEWMESHREG_DISCRETEGROUPCOFUNCTION_H

#include "DiscreteCostFunction.h"

namespace newmeshreg {

class DiscreteGroupCoCostFunction : public NonLinearSRegDiscreteCostFunction {

public:
    //---INITIALISATION---//
    void initialize(int numNodes, int numLabels, int numPairs, int numTriplets) override {
        m_num_nodes = numNodes;
        m_num_labels = numLabels;
        m_num_pairs = numPairs;
        m_num_triplets = numTriplets;
    }

    void set_warps(const std::vector<std::vector<newresampler::Mesh>>& ws) override {
        warps = ws;
    }

    void set_meshes(const std::vector<newresampler::Mesh>& source, const newresampler::Mesh& GRID, int num) override {
        _ORIG_MESHES = source;
        VERTICES_PER_SUBJ = GRID.nvertices();
        TRIPLETS_PER_SUBJ = GRID.ntriangles();
        _CONTROLMESHES.resize(num,GRID);
    }

    void reset_CPgrid(const newresampler::Mesh& grid, int num) override { _CONTROLMESHES[num] = grid; }
    void set_patch_data(const std::vector<std::map<int,float>>& patches) override { patch_data = patches; }
    void set_warp_rotations(const std::vector<std::vector<std::vector<NEWMAT::Matrix>>>& w_rot) override {
        warp_rotations = w_rot;
    }

    //---Compute costs---//
    double computeTripletCost(int triplet, int labelA, int labelB, int labelC) override;
    double computePairwiseCost(int pair, int labelA, int labelB) override;

private:
    std::vector<newresampler::Mesh> _CONTROLMESHES;
    std::vector<newresampler::Mesh> _ORIG_MESHES;
    std::vector<std::vector<newresampler::Mesh>> warps;

    std::vector<std::map<int,float>> patch_data;
    std::vector<std::vector<std::vector<NEWMAT::Matrix>>> warp_rotations;

    int TRIPLETS_PER_SUBJ = 0;
    int VERTICES_PER_SUBJ = 0;
};

} //namespace newmeshreg

#endif //NEWMESHREG_DISCRETEGROUPCOFUNCTION_H
