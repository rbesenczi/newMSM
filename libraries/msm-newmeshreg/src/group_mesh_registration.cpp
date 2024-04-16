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
#include "group_mesh_registration.h"

namespace newmeshreg {

void Group_Mesh_registration::initialize_level(int current_lvl) {

    check();
    if(cost[current_lvl] == "RIGID" || cost[current_lvl] == "AFFINE")
        throw MeshregException("AFFINE/RIGID registration is not supported in groupwise mode.");

    const std::vector<double> sigma(num_subjects, _sigma_in[current_lvl]);

    FEAT = std::make_shared<featurespace>(DATAlist);
    FEAT->set_smoothing_parameters(sigma);
    FEAT->set_cutthreshold(_threshold);
    FEAT->varnorm(_varnorm);
    FEAT->intensitynormalise(_IN, _cut);
    FEAT->is_sparse(_issparse);
    FEAT->set_nthreads(_numthreads);
    SPH_orig = FEAT->initialise(_genesis[current_lvl], MESHES, _exclude);

    PARAMETERS.insert(parameterPair("multivariate", FEAT->get_dim() > 1));

    newresampler::Mesh control = newresampler::make_mesh_from_icosa(_gridres[current_lvl]);
    newresampler::recentre(control);
    newresampler::true_rescale(control, RAD);

    model = std::make_shared<DiscreteGroupModel>(PARAMETERS);
    if(_debug) model->set_debug();
    model->set_featurespace(FEAT);
    model->set_meshspace(target_space, SPH_orig, num_subjects);
    if (is_masked) model->set_masks(mask);
    model->Initialize(control);
}

void Group_Mesh_registration::evaluate() {

    if (level == 1)
        ALL_SPH_REG.resize(num_subjects, SPH_orig);
    else
        for (int subject = 0; subject < num_subjects; subject++)
            ALL_SPH_REG[subject] = project_CPgrid(SPH_orig, ALL_SPH_REG[subject], subject);

    run_discrete_opt();

    if(_verbose) std::cout << "Exit main algorithm at level " << level << '.' << std::endl;
}

void Group_Mesh_registration::run_discrete_opt() {

    double energy = 0.0, newenergy = 0.0;
    int max_iter = std::get<int>(PARAMETERS.find("iters")->second);

    std::vector<newresampler::Mesh> previous_controlgrids(num_subjects);

    for (int subject = 0; subject < num_subjects; subject++)
        previous_controlgrids[subject] = model->get_CPgrid(subject);

    for (int iter = 0; iter < max_iter; iter++)
    {
        model->setupCostFunctionWeighting(combine_weighting());
        model->setupCostFunction();

#ifdef HAS_HOCR
        newenergy = Fusion::optimize(model, _verbose, _numthreads);
#else
        throw MeshregException("Groupwise mode is only supported in the HOCR version of MSM.");
#endif

        if (iter > 1 && /*iter % 2 != 0 && */ (energy-newenergy < newenergy * 0.01))
        {
            if (_verbose)
                std::cout << iter+1 << " iters has converged.\n"
                          << "New energy==" << newenergy << "\tPrevious energy==" << energy
                          << "\tEnergy decrease==" << energy - newenergy << std::endl;
            break;
        }

        if (iter > 1 && _verbose)
            std::cout << "New energy==" << newenergy << "\tPrevious energy==" << energy
                      << "\tEnergy decrease==" << energy - newenergy << std::endl;

        model->applyLabeling();

        for(int subject = 0; subject < num_subjects; subject++)
        {
            newresampler::Mesh transformed_controlgrid = model->get_CPgrid(subject);
            unfold(transformed_controlgrid, _verbose);
            newresampler::sphere_project_warp(ALL_SPH_REG[subject], previous_controlgrids[subject], transformed_controlgrid, _numthreads);
            unfold(ALL_SPH_REG[subject], _verbose);
            previous_controlgrids[subject] = transformed_controlgrid;
            model->reset_CPgrid(transformed_controlgrid, subject);
            model->reset_meshspace(ALL_SPH_REG[subject], subject);
        }
        energy = newenergy;
    }
}

void Group_Mesh_registration::transform(const std::string &filename) {
    for(int subject = 0; subject < num_subjects; ++subject) {
        newresampler::sphere_project_warp(MESHES[subject], SPH_orig, ALL_SPH_REG[subject], _numthreads);
        MESHES[subject].save(filename + "sphere-" + std::to_string(subject) + ".reg" + _surfformat);
    }
}

void Group_Mesh_registration::save_transformed_data(const std::string &filename) {
    for(int subject = 0; subject < num_subjects; ++subject) {
        std::shared_ptr<MISCMATHS::BFMatrix> data;
        set_data(DATAlist[subject], data, MESHES[subject]);
        newresampler::metric_resample(MESHES[subject], target_space, _numthreads).save(filename + "transformed_and_reprojected-" + std::to_string(subject) + _dataformat);
    }
}

} //namespace newmeshreg
