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
#include "DiscreteGroupCostFunction.h"

namespace newmeshreg {

double DiscreteGroupCostFunction::computeTripletCost(int triplet, int labelA, int labelB, int labelC) {

    int subject = std::floor((double)triplet / TRIPLETS_PER_SUBJ);

    int vertex_0 = _triplets[3*triplet  ] - subject * VERTICES_PER_SUBJ;
    int vertex_1 = _triplets[3*triplet+1] - subject * VERTICES_PER_SUBJ;
    int vertex_2 = _triplets[3*triplet+2] - subject * VERTICES_PER_SUBJ;

    newresampler::Triangle TRI((*ROTATIONS)[_triplets[3*triplet  ]] * _labels[labelA],
                               (*ROTATIONS)[_triplets[3*triplet+1]] * _labels[labelB],
                               (*ROTATIONS)[_triplets[3*triplet+2]] * _labels[labelC],0);
    newresampler::Triangle TRI_noDEF(_CONTROLMESHES[subject].get_coord(vertex_0),
                                     _CONTROLMESHES[subject].get_coord(vertex_1),
                                     _CONTROLMESHES[subject].get_coord(vertex_2), 0);

    if ((TRI.normal() | TRI_noDEF.normal()) < 0.0) return FOLDING;  // folded

    newresampler::Triangle TRI_ORIG(_ORIG_MESHES[subject].get_coord(vertex_0),
                                    _ORIG_MESHES[subject].get_coord(vertex_1),
                                    _ORIG_MESHES[subject].get_coord(vertex_2), 0);

    double strain_energy = calculate_triangular_strain(TRI_ORIG,TRI,_mu,_kappa,
                                                       std::shared_ptr<NEWMAT::ColumnVector>(), _k_exp);

    if(fixnan && std::isnan(strain_energy)) return FIX_NAN;
    else return _reglambda * MISCMATHS::pow(strain_energy,_rexp);
}

double DiscreteGroupCostFunction::computePairwiseCost(int pair, int labelA, int labelB) {

    std::vector<std::vector<double>> patch_data_A, patch_data_B;
    std::vector<double> weights;
    patch_data_A.reserve(1500);
    patch_data_B.reserve(1500);
    weights.reserve(1500);
    double pair_cost = 0.0;

    int subject_A = std::floor((double)_pairs[2*pair  ] / VERTICES_PER_SUBJ);
    int subject_B = std::floor((double)_pairs[2*pair+1] / VERTICES_PER_SUBJ);

    int node_A = _pairs[2*pair  ] - subject_A * VERTICES_PER_SUBJ;
    int node_B = _pairs[2*pair+1] - subject_B * VERTICES_PER_SUBJ;

    const std::map<int,std::vector<double>>& patchA = patch_data[subject_A * VERTICES_PER_SUBJ * m_num_labels + node_A * m_num_labels + labelA];
    const std::map<int,std::vector<double>>& patchB = patch_data[subject_B * VERTICES_PER_SUBJ * m_num_labels + node_B * m_num_labels + labelB];

    for (const auto& e: patchA) {
        auto it = patchB.find(e.first);
        if (it != patchB.end()) {
            patch_data_A.push_back(e.second);
            patch_data_B.push_back(it->second);
            if (is_masked) weights.push_back(std::abs(_MASK.get_pvalue(e.first)));
            else weights.push_back(1.0);
        }
    }

    int patch_size = patch_data_A.size();
    int dimensions = patch_data_A[0].size();

    std::vector<double> data_a(patch_size), data_b(patch_size);
    for(int dim = 0; dim < dimensions; dim++) {
        for (int datapoint = 0; datapoint < patch_data_A.size(); datapoint++) {
            data_a[datapoint] = patch_data_A[datapoint][dim];
            data_b[datapoint] = patch_data_B[datapoint][dim];
        }
        pair_cost += sim.get_sim_for_min(data_a, data_b, weights);
    }

    pair_cost /= dimensions;

    if(fixnan && std::isnan(pair_cost)) return FIX_NAN;
    else return pair_cost;
}

} //namespace newmeshreg
