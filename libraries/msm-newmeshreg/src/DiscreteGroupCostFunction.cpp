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

    if((TRI.normal() | TRI_noDEF.normal()) < 0) { return 1e7; } // folded

    newresampler::Triangle TRI_ORIG(_ORIG_MESHES[subject].get_coord(vertex_0),
                                    _ORIG_MESHES[subject].get_coord(vertex_1),
                                    _ORIG_MESHES[subject].get_coord(vertex_2), 0);

    return _reglambda * MISCMATHS::pow(calculate_triangular_strain(TRI_ORIG,TRI,_mu,_kappa,
                                           std::shared_ptr<NEWMAT::ColumnVector>(), _k_exp),_rexp);
}

double DiscreteGroupCostFunction::computePairwiseCost(int pair, int labelA, int labelB) {

    std::vector<double> patch_data_A, patch_data_B;
    patch_data_A.reserve(1500);
    patch_data_B.reserve(1500);

    int subject_A = std::floor((double)_pairs[2*pair  ] / VERTICES_PER_SUBJ);
    int subject_B = std::floor((double)_pairs[2*pair+1] / VERTICES_PER_SUBJ);

    int node_A = _pairs[2*pair  ] - subject_A * VERTICES_PER_SUBJ;
    int node_B = _pairs[2*pair+1] - subject_B * VERTICES_PER_SUBJ;

    const std::map<int,double>& patchA = patch_data[subject_A * VERTICES_PER_SUBJ * m_num_labels + node_A * m_num_labels + labelA];
    const std::map<int,double>& patchB = patch_data[subject_B * VERTICES_PER_SUBJ * m_num_labels + node_B * m_num_labels + labelB];

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
