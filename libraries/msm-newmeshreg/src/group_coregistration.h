#ifndef NEWMESHREG_GROUPCOMESHREG_H
#define NEWMESHREG_GROUPCOMESHREG_H

#include "mesh_registration.h"
#include "DiscreteGroupCoModel.h"

namespace newmeshreg {

class Group_coregistration : public Mesh_registration {

    std::vector<newresampler::Mesh> PAIR_SPH_REG;
    std::vector<newresampler::Mesh> warpsA;
    std::vector<newresampler::Mesh> warpsB;
    int num_subjects_A = 0, num_subjects_B = 0;
    newresampler::Mesh templ;

public:

    void initialize_level(int current_lvl) override;
    void evaluate() override;
    void run_discrete_opt() override;
    void transform(const std::string& filename) override;
    void save_transformed_data(const std::string& filename) override;
    void init_warps(int level);

    inline void set_warps(const std::string& warps_A, const std::string& warps_B) {
        std::vector<std::string> meshlist_A = read_ascii_list(warps_A);
        num_subjects_A = meshlist_A.size();
        warpsA.clear();
        warpsA.resize(num_subjects_A);
        for (int subject = 0; subject < num_subjects_A; ++subject) {
            if(_verbose) std::cout << "warp_A #" << subject << " is " << meshlist_A[subject] << std::endl;
            warpsA[subject].load(meshlist_A[subject]);
            recentre(warpsA[subject]);
            true_rescale(warpsA[subject], RAD);
        }

        std::vector<std::string> meshlist_B = read_ascii_list(warps_B);
        num_subjects_B = meshlist_B.size();
        warpsB.clear();
        warpsB.resize(num_subjects_B);
        for (int subject = 0; subject < num_subjects_B; ++subject) {
            if(_verbose) std::cout << "warp_B #" << subject << " is " << meshlist_B[subject] << std::endl;
            warpsB[subject].load(meshlist_B[subject]);
            recentre(warpsB[subject]);
            true_rescale(warpsB[subject], RAD);
        }
    }

    inline void set_meshes(const std::string& mesh_A, const std::string& mesh_B) {
        MESHES.clear();
        MESHES.resize(2);
        MESHES[0].load(mesh_A);
        recentre(MESHES[0]);
        true_rescale(MESHES[0], RAD);
        MESHES[1].load(mesh_B);
        recentre(MESHES[1]);
        true_rescale(MESHES[1], RAD);
        if (_verbose)
            std::cout << "mesh_A is " << mesh_A << "\nmesh_B is " << mesh_B << std::endl;
    }

    inline void set_data(const std::string& mean_A, const std::string& mean_B) {
        DATAlist.clear();
        DATAlist.resize(2);
        DATAlist[0] = mean_A;
        DATAlist[1] = mean_B;
        if (_verbose)
            for (int subject = 0; subject < DATAlist.size(); ++subject)
                std::cout << "Data #" << subject << " is " << DATAlist[subject] << std::endl;
    }

    inline void set_template(const std::string& tmpl) {
        if(_verbose) std::cout << "Template is " << tmpl << std::endl;
        templ.load(tmpl);
        recentre(templ);
        true_rescale(templ,RAD);
    }

    inline void saveSPH_reg(const std::string& filename) const override {
        PAIR_SPH_REG[0].save(filename + "sphere-0.LR.reg" + _surfformat);
        PAIR_SPH_REG[1].save(filename + "sphere-1.LR.reg" + _surfformat);
    }
};

} //namespace newmeshreg

#endif //NEWMESHREG_GROUPCOMESHREG_H
