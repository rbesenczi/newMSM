#include "DiscreteGroupCoCostFunction.h"

namespace newmeshreg {

double DiscreteGroupCoCostFunction::computeTripletCost(int triplet, int labelA, int labelB, int labelC) {

    int group = std::floor((double)triplet / TRIPLETS_PER_SUBJ);

    int vertex_0 = _triplets[3*triplet  ] - group * VERTICES_PER_SUBJ;
    int vertex_1 = _triplets[3*triplet+1] - group * VERTICES_PER_SUBJ;
    int vertex_2 = _triplets[3*triplet+2] - group * VERTICES_PER_SUBJ;

    std::vector<double> subject_costs(warps[group].size());
    for (int subject = 0; subject < warps[group].size(); subject++)
    {
        newresampler::Triangle TRI(warp_rotations[group][subject][vertex_0] * _labels[labelA],
                                   warp_rotations[group][subject][vertex_1] * _labels[labelB],
                                   warp_rotations[group][subject][vertex_2] * _labels[labelC],0);
        newresampler::Triangle TRI_noDEF(warps[group][subject].get_coord(vertex_0),
                                         warps[group][subject].get_coord(vertex_1),
                                         warps[group][subject].get_coord(vertex_2), 0);

        if((TRI.normal() | TRI_noDEF.normal()) < 0) { return 1e7; } // folded

        newresampler::Triangle TRI_ORIG(_ORIG_MESHES[group].get_coord(vertex_0),
                                        _ORIG_MESHES[group].get_coord(vertex_1),
                                        _ORIG_MESHES[group].get_coord(vertex_2), 0);

        subject_costs[subject] = _reglambda * MISCMATHS::pow(calculate_triangular_strain(TRI_ORIG,TRI,_mu,_kappa,
                                                                std::shared_ptr<NEWMAT::ColumnVector>(), _k_exp),_rexp);
    }

    double mean_reg_cost = std::accumulate(subject_costs.begin(), subject_costs.end(), 0.0) / subject_costs.size();
    double max_reg_cost = *(std::max_element(subject_costs.begin(), subject_costs.end()));
/*
    std::cout <<
        "Mean regcost==" << mean_reg_cost <<
        " Max regcost==" << max_reg_cost << '\t';
    for(const auto& e: subject_costs)
        std::cout << e << ' ';
    std::cout << std::endl;
*/
    return _reglambda * mean_reg_cost + 0.2 * max_reg_cost;
}

double DiscreteGroupCoCostFunction::computePairwiseCost(int pair, int labelA, int labelB) {

    std::vector<double> patch_data_A, patch_data_B;
    patch_data_A.reserve(1500);
    patch_data_B.reserve(1500);

    int subject_A = std::floor((double)_pairs[2*pair  ] / VERTICES_PER_SUBJ);
    int subject_B = std::floor((double)_pairs[2*pair+1] / VERTICES_PER_SUBJ);

    int node_A = _pairs[2*pair  ] - subject_A * VERTICES_PER_SUBJ;
    int node_B = _pairs[2*pair+1] - subject_B * VERTICES_PER_SUBJ;

    const std::map<int,float>& patchA = patch_data[subject_A * VERTICES_PER_SUBJ * m_num_labels + node_A * m_num_labels + labelA];
    const std::map<int,float>& patchB = patch_data[subject_B * VERTICES_PER_SUBJ * m_num_labels + node_B * m_num_labels + labelB];

    for (const auto& e: patchA) {
        auto it = patchB.find(e.first);
        if (it != patchB.end()) {
            patch_data_A.push_back(e.second);
            patch_data_B.push_back(it->second);
        }
    }

    return sim.get_sim_for_min(patch_data_A, patch_data_B);
}

} //namespace newmeshreg
