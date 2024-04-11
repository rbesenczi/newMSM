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
#ifndef NEWMESHREG_GROUPMESHREG_H
#define NEWMESHREG_GROUPMESHREG_H

#include "mesh_registration.h"
#include "DiscreteGroupModel.h"
#include "newresampler/mesh.h"

namespace newmeshreg {

class Group_Mesh_registration : public Mesh_registration {

    std::vector<newresampler::Mesh> ALL_SPH_REG;
    newresampler::Mesh target_space;
    newresampler::Mesh mask;
    int num_subjects = 1;
    bool is_masked = false;

public:
    void initialize_level(int current_lvl) override;
    void evaluate() override;
    void run_discrete_opt() override;
    void transform(const std::string& filename) override;
    void save_transformed_data(const std::string& filename) override;

    inline void set_inputs(const std::string& s) {
        std::vector<std::string> meshlist = read_ascii_list(s);
        num_subjects = meshlist.size();
        MESHES.clear();
        MESHES.resize(num_subjects);
        for (int subject = 0; subject < num_subjects; ++subject) {
            if(_verbose) std::cout << "Mesh #" << subject << " is " << meshlist[subject] << std::endl;
            MESHES[subject].load(meshlist[subject]);
            recentre(MESHES[subject]);
            true_rescale(MESHES[subject], RAD);
        }
    }

    inline void set_mask(const std::string& mask_file){
        mask = target_space;
        mask.load(mask_file, false, false);
        is_masked = true;
    }

    inline void set_data_list(const std::string& s) {
        DATAlist = read_ascii_list(s);
        if (_verbose)
            for (int subject = 0; subject < DATAlist.size(); ++subject)
                std::cout << "Data #" << subject << " is " << DATAlist[subject] << std::endl;
    }

    inline void set_template(const std::string &M) {
        if(_verbose) std::cout << "Template is " << M << std::endl;
        target_space.load(M);
        recentre(target_space);
        true_rescale(target_space, RAD);
    }

    inline void saveSPH_reg(const std::string& filename) const override {
        for(int subject = 0; subject < num_subjects; ++subject)
            ALL_SPH_REG[subject].save(filename + "sphere-" + std::to_string(subject) + ".LR.reg" + _surfformat);
    }
};

} //namespace newmeshreg

#endif //NEWMESHREG_GROUPMESHREG_H
