#ifndef NEWMESHREG_DISCRETEGROUPCOMODEL_H
#define NEWMESHREG_DISCRETEGROUPCOMODEL_H

#include "DiscreteModel.h"
#include "DiscreteGroupCostFunction.h"
#include "DiscreteGroupCoCostFunction.h"

namespace newmeshreg {

class DiscreteGroupCoModel : public NonLinearSRegDiscreteModel {

    std::vector<newresampler::Mesh> m_datameshes;
    std::vector<newresampler::Mesh> m_controlmeshes;
    newresampler::Mesh m_template;
    std::vector<NEWMAT::ColumnVector> spacings;
    std::vector<newresampler::Mesh> warpsA;
    std::vector<newresampler::Mesh> warpsB;

    int m_num_subjects = 0;
    int control_grid_size = 0;
    int cp_triangles = 0;

public:
    DiscreteGroupCoModel() = default;
    explicit DiscreteGroupCoModel(myparam& p) {
        set_parameters(p);
        costfct = std::make_shared<DiscreteGroupCoCostFunction>();
        costfct->set_parameters(p);
    }

    void set_warps(const std::vector<newresampler::Mesh>& warps_A, const std::vector<newresampler::Mesh>& warps_B) override {
        warpsA = warps_A;
        warpsB = warps_B;
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

#endif // NEWMESHREG_DISCRETEGROUPCOMODEL_H