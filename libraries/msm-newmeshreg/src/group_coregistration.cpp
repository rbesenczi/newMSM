#include "group_coregistration.h"

namespace newmeshreg {

void Group_coregistration::initialize_level(int current_lvl) {

    check();
    if(cost[current_lvl] == "RIGID" || cost[current_lvl] == "AFFINE")
        throw MeshregException("AFFINE/RIGID registration is not supported in groupwise co-registration mode yet.");

    const std::vector<double> sigma(2, _sigma_in[current_lvl]);

    FEAT = std::make_shared<featurespace>(DATAlist);
    FEAT->set_smoothing_parameters(sigma);
    FEAT->set_cutthreshold(_threshold);
    FEAT->varnorm(_varnorm);
    FEAT->intensitynormalize(_IN, _cut);
    FEAT->is_sparse(_issparse);
    FEAT->set_nthreads(_numthreads);
    SPH_orig = FEAT->initialize(_genesis[current_lvl], MESHES, _exclude);
    if(FEAT->get_dim() > 1)
        throw MeshregException("Multivariate registration is not supported in groupwise co-registration mode yet.");
    PARAMETERS.insert(parameterPair("multivariate", false));

    newresampler::Mesh control = newresampler::make_mesh_from_icosa(_gridres[current_lvl]);
    newresampler::recentre(control);
    newresampler::true_rescale(control, RAD);

    control_warps = init_warps(current_lvl);

    model = std::make_shared<DiscreteGroupCoModel>(PARAMETERS);
    if(_debug) model->set_debug();
    model->set_featurespace(FEAT);
    model->set_meshspace(templ, SPH_orig, 2);
    model->Initialize(control);
}

void Group_coregistration::evaluate() {

    if (level == 1)
        PAIR_SPH_REG.resize(2, SPH_orig);
    else
        for (int subject = 0; subject < 2; subject++) {
            PAIR_SPH_REG[subject] = project_CPgrid(SPH_orig, PAIR_SPH_REG[subject], subject);
            for (int warp = 0; warp < warps[subject].size(); warp++)
                control_warps[subject][warp] = project_CPgrid(SPH_orig, control_warps[subject][warp], subject);
        }

    run_discrete_opt();

    if(_verbose) std::cout << "Exit main algorithm at level " << level << '.' << std::endl;
}

void Group_coregistration::run_discrete_opt() {

    double energy = 0.0, newenergy = 0.0;

    std::vector<newresampler::Mesh> previous_controlgrids(2);

    for (int subject = 0; subject < 2; subject++)
        previous_controlgrids[subject] = model->get_CPgrid(subject);

    for(int iter = 1; iter <= std::get<int>(PARAMETERS.find("iters")->second); iter++)
    {
        model->set_warps(control_warps);
        model->setupCostFunction();

#ifdef HAS_HOCR
        newenergy = Fusion::optimize(model, _verbose, _numthreads);
#else
        throw MeshregException("Groupwise mode is only supported in the HOCR version of MSM.");
#endif

        if(iter > 1 && iter % 2 != 0 && (energy-newenergy < newenergy*0.01)) // convergence when less than 1% decrease in energy
        {
            if (_verbose)
                std::cout << iter << " level has converged.\n"
                          << "New energy==" << newenergy << "\tPrevious energy==" << energy
                          << "\tEnergy decrease==" << energy - newenergy << std::endl;
            break;
        }

        if (iter > 1 && _verbose) std::cout << "New energy==" << newenergy << "\tPrevious energy==" << energy
                                            << "\tEnergy decrease==" << energy - newenergy << std::endl;

        model->applyLabeling();

        for(int group = 0; group < 2; group++)
        {
            newresampler::Mesh transformed_controlgrid = model->get_CPgrid(group);
            unfold(transformed_controlgrid, _verbose);
            newresampler::barycentric_mesh_interpolation(PAIR_SPH_REG[group], previous_controlgrids[group], transformed_controlgrid, _numthreads);
            unfold(PAIR_SPH_REG[group], _verbose);

            for(int warp = 0; warp < warps[group].size(); ++warp)
            {
                newresampler::barycentric_mesh_interpolation(control_warps[group][warp], previous_controlgrids[group], transformed_controlgrid);
                unfold(control_warps[group][warp]);
            }

            previous_controlgrids[group] = transformed_controlgrid;
            model->reset_CPgrid(transformed_controlgrid, group);
            model->reset_meshspace(PAIR_SPH_REG[group], group);

        }
        energy = newenergy;
    }
}

void Group_coregistration::transform(const std::string &filename) {
    for(int subject = 0; subject < 2; ++subject)
    {
        newresampler::barycentric_mesh_interpolation(MESHES[subject], SPH_orig, PAIR_SPH_REG[subject], _numthreads);
        MESHES[subject].save(filename + "sphere-" + std::to_string(subject) + ".reg" + _surfformat);
    }
}

void Group_coregistration::save_transformed_data(const std::string &filename) {
    for(int subject = 0; subject < 2; ++subject)
    {
        std::shared_ptr<MISCMATHS::BFMatrix> data;
        newmeshreg::set_data(DATAlist[subject], data, MESHES[subject]);
        newresampler::metric_resample(MESHES[subject], templ, _numthreads).save(filename + "transformed_and_reprojected-" + std::to_string(subject) + _dataformat);
    }
    save_warps();
}

std::vector<std::vector<newresampler::Mesh>> Group_coregistration::init_warps(int level) {

    std::vector<std::vector<newresampler::Mesh>> control_warps(2);
    control_warps[0].resize(warps[0].size());
    control_warps[1].resize(warps[1].size());

    newresampler::Mesh new_ico = newresampler::make_mesh_from_icosa(_gridres[level]);
    newresampler::recentre(new_ico);
    newresampler::true_rescale(new_ico, RAD);

    for (int group = 0; group < 2; group++)
        for (int warp = 0; warp < warps[group].size(); warp++)
            control_warps[group][warp] = newresampler::surface_resample(warps[group][warp], MESHES[group], new_ico);

    return control_warps;
}

void Group_coregistration::save_warps() {

    for(int group = 0; group < 2; group++)
        for(int warp = 0; warp < warps[group].size(); warp++) {
            newresampler::barycentric_mesh_interpolation(warps[group][warp], SPH_orig, control_warps[group][warp], _numthreads);
            warps[group][warp].save(
                    _outdir + "sphere-" + std::to_string(group) + '.' + std::to_string(warp) + ".reg." + _surfformat);
        }
}

} //namespace newmeshreg
