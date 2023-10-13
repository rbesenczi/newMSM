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
        throw MeshregException("AFFINE/RIGID registration is not supported in groupwise mode yet.");

    const std::vector<double> sigma(num_subjects, _sigma_in[current_lvl]);

    FEAT = std::make_shared<featurespace>(DATAlist);
    FEAT->set_smoothing_parameters(sigma);
    FEAT->set_cutthreshold(_threshold);
    FEAT->varnorm(_varnorm);
    FEAT->intensitynormalize(_IN, _cut);
    FEAT->is_sparse(_issparse);
    FEAT->set_nthreads(_numthreads);
    SPH_orig = FEAT->initialize(_genesis[current_lvl], MESHES, _exclude);
    if(FEAT->get_dim() > 1)
        throw MeshregException("Multivariate registration is not supported in groupwise mode yet.");
    PARAMETERS.insert(parameterPair("multivariate", false));

    newresampler::Mesh control = newresampler::make_mesh_from_icosa(_gridres[current_lvl]);
    newresampler::recentre(control);
    newresampler::true_rescale(control, RAD);

    model = std::make_shared<DiscreteGroupModel>(PARAMETERS);
    if(_debug) model->set_debug();
    model->set_featurespace(FEAT);
    model->set_meshspace(templ, SPH_orig, num_subjects);
    model->Initialize(control);
}

void Group_Mesh_registration::evaluate() {

    if (level == 1)
        ALL_SPH_REG.resize(num_subjects, SPH_orig);
    else
        for (int subject = 0; subject < num_subjects; subject++)
            ALL_SPH_REG[subject] = project_CPgrid(SPH_orig, ALL_SPH_REG[subject], subject);

    run_discrete_opt(ALL_SPH_REG);

    if(_verbose) std::cout << "Exit main algorithm at level " << level << '.' << std::endl;
}

void Group_Mesh_registration::run_discrete_opt(std::vector<newresampler::Mesh>& meshes) {

    double energy = 0.0, newenergy = 0.0;
    std::vector<newresampler::Mesh> previous_controlgrids(num_subjects);

    for(int iter = 1; iter <= boost::get<int>(PARAMETERS.find("iters")->second); iter++) {

        for (int subject = 0; subject < num_subjects; subject++) {
            model->reset_meshspace(ALL_SPH_REG[subject], subject);
            previous_controlgrids[subject] = model->get_CPgrid(subject);
        }

        model->setupCostFunction();

        if(_debug)
            for(int subject = 0; subject < num_subjects; subject++) {
                newresampler::Mesh tmp = MESHES[subject];
                newresampler::barycentric_mesh_interpolation(tmp, SPH_orig, ALL_SPH_REG[subject], _numthreads);
                std::shared_ptr<MISCMATHS::BFMatrix> data;
                set_data(DATAlist[subject], data, tmp);
                newresampler::metric_resample(tmp, templ, _numthreads).save(
                        _outdir + "transformed_and_reprojected-" + std::to_string(subject)
                        + "before-level-" + std::to_string(level) + "-iter-" + std::to_string(iter) + _dataformat);
            }

#ifdef HAS_HOCR
        newenergy = Fusion::optimize(model, _verbose, _numthreads);
#else
        throw MeshregException("Groupwise mode is only supported in the HOCR version of MSM.");
#endif
        if(iter > 1 && iter % 2 != 0 && energy-newenergy < 1) {
            if (_verbose)
                std::cout << iter << " level has converged.\n"
                          << "New energy==" << newenergy << "\tPrevious energy==" << energy
                          << "\tEnergy decrease==" << energy - newenergy << std::endl;
            break;
        }

        if (iter > 1 && _verbose) std::cout << "\tEnergy decrease==" << energy - newenergy << std::endl;

        model->applyLabeling();

        for(int subject = 0; subject < num_subjects; subject++)
        {
            auto transformed_controlgrid = model->get_CPgrid(subject);
            unfold(transformed_controlgrid, _verbose);
            model->reset_CPgrid(transformed_controlgrid, subject);
            newresampler::barycentric_mesh_interpolation(ALL_SPH_REG[subject], previous_controlgrids[subject], transformed_controlgrid, _numthreads);
            unfold(ALL_SPH_REG[subject], _verbose);
            if(_debug) {
                newresampler::Mesh tmp = MESHES[subject];
                newresampler::barycentric_mesh_interpolation(tmp, SPH_orig, ALL_SPH_REG[subject], _numthreads);
                std::shared_ptr<MISCMATHS::BFMatrix> data;
                set_data(DATAlist[subject], data, tmp);
                newresampler::metric_resample(tmp, templ, _numthreads).save(
                        _outdir + "transformed_and_reprojected-" + std::to_string(subject)
                        + "after-level-" + std::to_string(level) + "-iter-" + std::to_string(iter) + _dataformat);
            }
        }

        energy = newenergy;
    }
}

void Group_Mesh_registration::transform(const std::string &filename) {
    for(int subject = 0; subject < num_subjects; ++subject)
    {
        newresampler::barycentric_mesh_interpolation(MESHES[subject], SPH_orig, ALL_SPH_REG[subject], _numthreads);
        MESHES[subject].save(filename + "sphere-" + std::to_string(subject) + ".reg" + _surfformat);
    }
}

void Group_Mesh_registration::save_transformed_data(const std::string &filename) {
    for(int subject = 0; subject < num_subjects; ++subject)
    {
        std::shared_ptr<MISCMATHS::BFMatrix> data;
        set_data(DATAlist[subject], data, MESHES[subject]);
        newresampler::metric_resample(MESHES[subject], templ, _numthreads).save(filename + "transformed_and_reprojected-" + std::to_string(subject) + _dataformat);
    }
}

} //namespace newmeshreg
