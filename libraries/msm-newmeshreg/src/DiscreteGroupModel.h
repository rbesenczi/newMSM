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
#ifndef NEWMESHREG_DISCRETEGROUPMODEL_H
#define NEWMESHREG_DISCRETEGROUPMODEL_H

#include "DiscreteModel.h"
#include "DiscreteGroupCostFunction.h"

namespace newmeshreg {

class DiscreteGroupModel : public NonLinearSRegDiscreteModel {

    std::vector<newresampler::Mesh> m_datameshes;
    std::vector<newresampler::Mesh> m_controlmeshes;
    newresampler::Mesh m_template;
    std::vector<std::map<int,double>> patch_data;
    std::vector<NEWMAT::ColumnVector> spacings;

    int m_num_subjects = 0;
    int control_grid_size = 0;
    int cp_triangles = 0;

public:
    DiscreteGroupModel() = default;
    explicit DiscreteGroupModel(myparam& p) {
        set_parameters(p);
        costfct = std::make_shared<DiscreteGroupCostFunction>();
        costfct->set_parameters(p);
    }

    void set_meshspace(const newresampler::Mesh& target, const newresampler::Mesh& source, int num) override {
        m_template = target;
        m_datameshes.clear();
        m_datameshes.resize(num, source);
        m_num_subjects = num;
    }

    void reset_meshspace(const newresampler::Mesh& source, int num) override {
        m_datameshes[num] = source;
    }

    void reset_CPgrid(const newresampler::Mesh& grid, int num) override {
        m_controlmeshes[num] = grid;
        costfct->reset_CPgrid(grid, num);
    }

    void warp_CPgrid(newresampler::Mesh& start, newresampler::Mesh& end, int num) override {
        newresampler::barycentric_mesh_interpolation(m_controlmeshes[num], start, end, _nthreads);
        unfold(m_controlmeshes[num], m_verbosity);
    }

    void applyLabeling() override {
        for (int subject = 0; subject < m_num_subjects; subject++)
            for (int vertex = 0; vertex < control_grid_size; vertex++)
                m_controlmeshes[subject].set_coord(vertex, m_ROT[subject * control_grid_size + vertex] * m_labels[labeling[subject * control_grid_size + vertex]]);
    }

    newresampler::Mesh get_CPgrid(int num) override { return m_controlmeshes[num]; }

    void get_spacings();
    void initialize_pairs();
    void estimate_pairs() override;
    void estimate_triplets() override;
    void get_patch_data();

    void Initialize(const newresampler::Mesh& controlgrid) override;
    void get_rotations() override;
    void setupCostFunction() override;
};

} //namespace newmeshreg

#endif //NEWMESHREG_DISCRETEGROUPMODEL_H
