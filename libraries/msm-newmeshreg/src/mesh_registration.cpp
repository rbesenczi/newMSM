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
#include "mesh_registration.h"

namespace newmeshreg {

Mesh_registration::Mesh_registration(){
    MESHES.resize(2, newresampler::Mesh());
}

void Mesh_registration::run_multiresolutions(const std::string& parameters) {

    parse_reg_options(parameters);

    if(_verbose) std::cout << "Starting multiresolution with " << _resolutionlevels << " levels." << std::endl;

    for(int i = 0; i < _resolutionlevels; i++)
    {
        level = i + 1;
        if(_verbose) std::cout << "Initialising level " << level << std::endl;

        fix_parameters_for_level(i);
        initialize_level(i);
        evaluate();
        if(isrigid) rigidcf.reset();
    }

    transform(_outdir);
    saveSPH_reg(_outdir);
    save_transformed_data(_outdir);
}

void Mesh_registration::initialize_level(int current_lvl) {

    check();

    FEAT = std::make_shared<featurespace>(CMfile_in, CMfile_ref);
    FEAT->set_smoothing_parameters({_sigma_in[current_lvl], _sigma_ref[current_lvl]});
    FEAT->set_cutthreshold(_threshold); // will also generate exclusion masks at the same mesh resolution as datagrid
    FEAT->varnorm(_varnorm);// variance normalises
    FEAT->intensitynormalize(_IN, _cut); // matches the intensities of the source to the target (will rescale all to the top feature of the target if scale is true)
    FEAT->is_sparse(_issparse);
    FEAT->set_nthreads(_numthreads);
    SPH_orig = FEAT->initialize(_genesis[current_lvl], MESHES, _exclude);  // downsamples and smooths data, creates and exclusion mask if exclude is true

    SPHin_CFWEIGHTING = downsample_cfweighting(_sigma_in[current_lvl], SPH_orig, IN_CFWEIGHTING, FEAT->get_input_excl());
    SPHref_CFWEIGHTING = downsample_cfweighting(_sigma_ref[current_lvl], SPH_orig, REF_CFWEIGHTING, FEAT->get_reference_excl());

    if(cost[current_lvl] == "RIGID" || cost[current_lvl] == "AFFINE")
    {
        rigidcf = std::make_shared<Rigid_cost_function>(SPH_orig, SPH_orig, FEAT);
        rigidcf->set_parameters(PARAMETERS);
        if(_simval[current_lvl] != 1) rigidcf->set_simmeasure(1);
        rigidcf->Initialize();
        isrigid = true;
    }

    if(cost[current_lvl] == "DISCRETE")
    {
        isrigid = false;

        PARAMETERS.insert(parameterPair("multivariate",FEAT->get_dim() > 1));

        model = std::make_shared<NonLinearSRegDiscreteModel>(PARAMETERS);

        if(_debug) model->set_debug();
        model->set_featurespace(FEAT);
        model->set_meshspace(SPH_orig, SPH_orig, 0);
        newresampler::Mesh CONTROL = newresampler::make_mesh_from_icosa(_gridres[current_lvl]);
        newresampler::recentre(CONTROL);
        newresampler::true_rescale(CONTROL, RAD);

        if(_anat)
        {
            std::vector<std::vector<int>> ANAT_face_neighbourhood;
            std::vector<std::map<int,double>> ANAT_to_CP_baryweights;
            newresampler::Mesh aICO = resample_anatomy(CONTROL, ANAT_to_CP_baryweights, ANAT_face_neighbourhood, current_lvl);
            newresampler::Mesh ANAT_target = newresampler::surface_resample(ref_anat,MESHES[1],aICO, _numthreads);
            model->set_anatomical_meshspace(aICO, ANAT_target, aICO, ANAT_orig);
            model->set_anatomical_neighbourhood(ANAT_to_CP_baryweights, ANAT_face_neighbourhood);
        }
        else
        {
            if (_regmode == 5) throw MeshregException("STRAINS based regularisation requires anatomical meshes");
            else if (_regmode == 4) throw MeshregException("You have specified angular based penalisation of anatomical warp, for which you require anatomical meshes");
        }

        model->Initialize(CONTROL);
    }
}

newresampler::Mesh Mesh_registration::resample_anatomy(const newresampler::Mesh& control_grid,
                                                    std::vector<std::map<int,double>>& baryweights,
                                                    std::vector<std::vector<int>>& ANAT_to_CPgrid_neighbours,
                                                    int current_lvl) {

    if (in_anat.nvertices() != MESHES[0].nvertices() || ref_anat.nvertices() != MESHES[1].nvertices())
        throw MeshregException("MeshREG ERROR:: input/reference anatomical mesh resolution is inconsistent with input/reference spherical mesh resolution.");

    newresampler::Mesh ANAT_ico = control_grid;

    // save face neighbourhood relationships during retesselation
    std::vector<std::vector<int>> FACE_neighbours, FACE_neighbours_tmp;
    std::vector<int> tmp;

    if((_anatres[current_lvl] - _gridres[current_lvl]) > 0)
    {
        for(int i = 0; i < (_anatres[current_lvl] - _gridres[current_lvl]); i++)
        {
            newresampler::retessellate(ANAT_ico, FACE_neighbours_tmp);
            if(i > 0)
            {
                ANAT_to_CPgrid_neighbours.clear();
                for(unsigned int j=0;j<FACE_neighbours.size();j++)
                {
                    ANAT_to_CPgrid_neighbours.push_back(tmp);
                    for(unsigned int k=0;k<FACE_neighbours[j].size();k++)
                    {
                        ANAT_to_CPgrid_neighbours[j].insert(ANAT_to_CPgrid_neighbours[j].begin(),
                                FACE_neighbours_tmp[FACE_neighbours[j][k]].begin(),
                                FACE_neighbours_tmp[FACE_neighbours[j][k]].end());
                    }
                }
                FACE_neighbours=ANAT_to_CPgrid_neighbours;
            }
            else
            {
                FACE_neighbours=FACE_neighbours_tmp;
            }
            FACE_neighbours_tmp.clear();
        }

        if (ANAT_to_CPgrid_neighbours.empty())
            ANAT_to_CPgrid_neighbours = FACE_neighbours;
            //i.e. for one increase in resolution use result directly from restesselation
    }
    else
    {
        for (int i = 0; i < control_grid.ntriangles(); i++)
        {
            ANAT_to_CPgrid_neighbours.push_back(tmp);
            ANAT_to_CPgrid_neighbours[i].push_back(i);
        }
    }

    true_rescale(ANAT_ico,RAD);

    // now get barycentric weights
    baryweights.resize(ANAT_ico.nvertices(),std::map<int,double>());

    for(int i = 0; i < ANAT_to_CPgrid_neighbours.size(); i++)
    {
        int id0 = control_grid.get_triangle(i).get_vertex_no(0);
        int id1 = control_grid.get_triangle(i).get_vertex_no(1);
        int id2 = control_grid.get_triangle(i).get_vertex_no(2);
        newresampler::Point v0 = control_grid.get_triangle_vertex(i,0);
        newresampler::Point v1 = control_grid.get_triangle_vertex(i,1);
        newresampler::Point v2 = control_grid.get_triangle_vertex(i,2);

        for(const int& j : ANAT_to_CPgrid_neighbours[i])
        {
            for(int k = 0; k < 3; k++)
            {
                const newresampler::Point& ci = ANAT_ico.get_triangle_vertex(j,k);
                int id = ANAT_ico.get_triangle(j).get_vertex_no(k);
                baryweights[id] = calc_barycentric_weights(v0, v1, v2, ci, id0, id1, id2);
            }
        }
    }

    ANAT_orig = newresampler::surface_resample(in_anat,MESHES[0],ANAT_ico, _numthreads);

    return ANAT_ico;
}

NEWMAT::Matrix Mesh_registration::downsample_cfweighting(double sigma,
                                                         const newresampler::Mesh& SPH,
                                                         std::shared_ptr<newresampler::Mesh> CFWEIGHTING,
                                                         std::shared_ptr<newresampler::Mesh> EXCL) {
    NEWMAT::Matrix newdata;

    if(EXCL)
    {
        if (!CFWEIGHTING)
            CFWEIGHTING = std::make_shared<newresampler::Mesh>(*EXCL);

        newdata = newresampler::nearest_neighbour_interpolation(*CFWEIGHTING, SPH, _numthreads, EXCL).get_pvalues();
    }
    else if(CFWEIGHTING)
    {
        newdata = newresampler::nearest_neighbour_interpolation(*CFWEIGHTING, SPH, _numthreads, EXCL).get_pvalues();
    }
    else
    {
        newdata.ReSize(1, SPH.nvertices());
        newdata = 1;
    }

    return newdata;
}

//---MAIN FUNCTION---//
void Mesh_registration::evaluate() {
    //Initialise deformation mesh
    SPH_reg = project_CPgrid(SPH_orig,SPH_reg);
    // first project data grid through any predefined transformation or,
    // from transformation from previous resolution level.

    if(isrigid)
    {
        rigidcf->update_source(SPH_reg);
        SPH_reg = rigidcf->run();
    }
    else
    {
        run_discrete_opt();
    }

    if(_verbose) std::cout << "Exit main algorithm." << std::endl;
}

void Mesh_registration::transform(const std::string& filename) {

    barycentric_mesh_interpolation(MESHES[0], SPH_orig, SPH_reg, _numthreads);
    MESHES[0].save(filename + "sphere.reg" + _surfformat);
}

void Mesh_registration::save_transformed_data(const std::string& filename) {

    if(_verbose) std::cout << "Saving and transforming data." << std::endl;

    std::string path = filename + "transformed_and_reprojected" + _dataformat;

    std::shared_ptr<MISCMATHS::BFMatrix> DATA, DATAREF;
    std::shared_ptr<newresampler::Mesh> IN_EXCL, REF_EXCL;

    // binarize costfunction weights for use as exclusion masks during resampling
    set_data(CMfile_in,DATA,MESHES[0]);
    set_data(CMfile_ref,DATAREF,MESHES[1]);

    if(_exclude)
    { // no longer necessary as EXCL masks aren't downsampled anymore? Could save from initialisation
        IN_EXCL= std::make_shared<newresampler::Mesh>(create_exclusion(MESHES[0],_threshold[0],_threshold[1]));
        REF_EXCL= std::make_shared<newresampler::Mesh>(create_exclusion(MESHES[1],_threshold[0],_threshold[1]));
    }
    if(_IN)
    { // INTENSITY NORMALIZE
        if(_verbose) std::cout << "Intensity normalise." << std::endl;
        multivariate_histogram_normalization(*DATA,*DATAREF,IN_EXCL,REF_EXCL, _numthreads);
    }

    MESHES[0].set_pvalues(DATA->AsMatrix());
    newresampler::Mesh tmp;

    tmp = newresampler::metric_resample(MESHES[0], MESHES[1], _numthreads, IN_EXCL);

    DATA = std::make_shared<MISCMATHS::FullBFMatrix>(tmp.get_pvalues());

    newresampler::Mesh TRANSFORMED = MESHES[1];
    std::shared_ptr<MISCMATHS::FullBFMatrix > pin = std::dynamic_pointer_cast<MISCMATHS::FullBFMatrix>(DATA);

    if(pin)
    {
      TRANSFORMED.set_pvalues(DATA->AsMatrix());
      TRANSFORMED.save(path);
    }
    else
        DATA->Print(path);

    if(_anat)
    {
        newresampler::Mesh ANAT_TRANS = newresampler::project_mesh(MESHES[0], MESHES[1],ref_anat, _numthreads);
        ANAT_TRANS.save(_outdir + "anat.reg.surf");

        in_anat.estimate_normals();
        ANAT_TRANS.estimate_normals();
        if(_verbose) std::cout << "Calculate strains." << std::endl;
        newresampler::Mesh STRAINSmesh = calculate_strains(2, in_anat, ANAT_TRANS, _numthreads);
        STRAINSmesh.save(_outdir + "STRAINS.func");
    }
}

//---PROJECT TRANSFORMATION FROM PREVIOUS LEVEL TO UPSAMPLED SOURCE---//
newresampler::Mesh Mesh_registration::project_CPgrid(newresampler::Mesh SPH_in, const newresampler::Mesh& REG, int num) {
    // num indices which warp for group registration

    if(level == 1)
    {
        if(transformed_mesh.nvertices() > 0)
        { // project into alignment with transformed mesh
            if (transformed_mesh == MESHES[num]) std::cout << " WARNING:: transformed mesh has the same coordinates as the input mesh " << std::endl;
            else
            {
                barycentric_mesh_interpolation(SPH_in, MESHES[num], transformed_mesh, _numthreads);
                if (model) model->warp_CPgrid(MESHES[num],transformed_mesh, num);
                // for tri clique model control grid is continously deformed
            }
        }
    }
    else
    {   // following first round always start by projecting Control and data grids through warp at previous level
        // PROJECT CPgrid into alignment with warp from previous level
        newresampler::Mesh icotmp = newresampler::make_mesh_from_icosa(REG.get_resolution());
        true_rescale(icotmp,RAD);
        // project datagrid though warp defined for the high resolution meshes (the equivalent to if registration is run one level at a time )
        newresampler::Mesh inorig = MESHES[num], incurrent = MESHES[num];
        barycentric_mesh_interpolation(incurrent,icotmp,REG, _numthreads);
        barycentric_mesh_interpolation(SPH_in,inorig,incurrent, _numthreads);
        if(model) model->warp_CPgrid(inorig, incurrent, num);
        if(_debug) incurrent.save(_outdir + "sphere.regLR.Res" + std::to_string(level) + ".surf");
    }

    unfold(SPH_in, _verbose);

    return SPH_in;
}

void Mesh_registration::run_discrete_opt() {

    double energy = 0.0, newenergy = 0.0;

    for(int iter = 1; iter <= std::get<int>(PARAMETERS.find("iters")->second); iter++) {
        // resample and combine the reference cost function weighting with the source if provided
        NEWMAT::Matrix CombinedWeight;
        if(_incfw && _refcfw) {
            NEWMAT::Matrix ResampledRefWeight = SPHref_CFWEIGHTING;
            newresampler::Mesh targetmesh = model->get_TARGET();
            targetmesh.set_pvalues(ResampledRefWeight);
            ResampledRefWeight = newresampler::metric_resample(targetmesh, SPH_reg, _numthreads).get_pvalues();
            CombinedWeight = combine_costfunction_weighting(SPHin_CFWEIGHTING, ResampledRefWeight);
        }
        else {
            CombinedWeight.resize(1, SPH_reg.nvertices());
            CombinedWeight = 1;
        }
        model->setupCostFunctionWeighting(CombinedWeight);

        model->reset_meshspace(SPH_reg,0);
        model->setupCostFunction();

        if(_verbose) std::cout << "Run optimisation." << std::endl;

        if(_discreteOPT == "MCMC") {
            if(!_tricliquelikeihood) model->computeUnaryCosts();
            newenergy = MCMC::optimise(model, _verbose, _mciters[level-1]);
        }
        else if(_discreteOPT == "FastPD") {
#ifdef HAS_FPD
            model->computeUnaryCosts();
            model->computePairwiseCosts();
            FPD::FastPD opt(model, 100);
            newenergy = opt.run();
            opt.getLabeling(model->getLabeling());
#else
            throw MeshregException("FastPD is not supported in this version of newMSM. Please use MCMC.");
#endif
        }
        else if(_discreteOPT == "HOCR") {
#ifdef HAS_HOCR
            newenergy = Fusion::optimize(model, _verbose, _numthreads);
#else
#ifdef HAS_FPD
            throw MeshregException("HOCR is not supported in this version of newMSM. Please use MCMC or FastPD.");
#else
            throw MeshregException("HOCR and FastPD are not supported in this version of newMSM. Please use MCMC.");
#endif
#endif
        }
        else
            throw MeshregException("Unrecognized optimiser");

        if(iter > 1 && ((iter - 1) % 2 == 0) && (energy - newenergy < 0.001) && _discreteOPT != "MCMC") {
            if(_verbose) {
                std::cout << iter << " level has converged." << std::endl;
                std::cout <<  "newenergy " << newenergy <<  "\tenergy " << energy
                          <<  "\tEnergy decrease: " <<  energy-newenergy << std::endl;
            }
            break;
        }

        if(_verbose) {
            std::cout << "newenergy " << newenergy << "\tenergy " << energy
                      << "\tEnergy decrease: " << energy - newenergy << std::endl;
        }

        newresampler::Mesh previous_controlgrid = model->get_CPgrid(0);

        model->applyLabeling();
        // apply these choices in order to deform the CP grid
        newresampler::Mesh transformed_controlgrid = model->get_CPgrid(0);
        // use the control point updates to warp the source mesh
        newresampler::barycentric_mesh_interpolation(SPH_reg, previous_controlgrid, transformed_controlgrid, _numthreads);
        // higher order frameowrk continuous deforms the CP grid whereas the original FW resets the grid each time
        unfold(transformed_controlgrid, _verbose);
        model->reset_CPgrid(transformed_controlgrid,0); // source mesh is updated and control point grids are reset
        unfold(SPH_reg, _verbose);
        energy = newenergy;
    }
}

void Mesh_registration::set_input(const newresampler::Mesh& M) {
    MESHES[0] = M;
    recentre(MESHES[0]);
    true_rescale(MESHES[0], RAD);
}

void Mesh_registration::set_input(const std::string &M) {
    MESHES[0].load(M);
    recentre(MESHES[0]);
    true_rescale(MESHES[0], RAD);
}

void Mesh_registration::set_reference(const newresampler::Mesh &M) {
    MESHES[1] = M;
    recentre(MESHES[1]);
    true_rescale(MESHES[1], RAD);
}

void Mesh_registration::set_reference(const std::string& M) {
    MESHES[1].load(M);
    recentre(MESHES[1]);
    true_rescale(MESHES[1], RAD);
}

void Mesh_registration::set_anatomical(const std::string &M1, const std::string &M2) {
    _anat = true;
    in_anat.load(M1);
    ref_anat.load(M2);
}

void Mesh_registration::set_transformed(const std::string &M) {
    transformed_mesh.load(M);
    true_rescale(MESHES[1], RAD);
}

void Mesh_registration::set_input_cfweighting(const std::string& E) {
    IN_CFWEIGHTING = std::make_shared<newresampler::Mesh>(MESHES[0]);
    IN_CFWEIGHTING->load(E, false, false);
    true_rescale(*IN_CFWEIGHTING, RAD);
    _incfw = true;
}

void Mesh_registration::set_reference_cfweighting(const std::string& E) {
    REF_CFWEIGHTING = std::make_shared<newresampler::Mesh>(MESHES[1]);
    REF_CFWEIGHTING->load(E, false, false);
    true_rescale(*REF_CFWEIGHTING, RAD);
    _refcfw = true;
}

void Mesh_registration::parse_reg_options(const std::string &parameters)
{
    std::string title = "newmsm configuration parameters";
    std::string examples;
    Utilities::OptionParser options(title,examples);

    std::vector<std::string> costdefault;
    Utilities::Option<std::vector<std::string>> optimizer(std::string("--opt"),costdefault,
                                               std::string("optimisation method. Choice of: AFFINE/RIGID, DISCRETE (default)"),
                                     false,Utilities::requires_argument);
    std::vector<int> intdefault;
    Utilities::Option<std::vector<int>> simval(std::string("--simval"), intdefault,
                                    std::string("code for determining which similarty measure is used to assess cost during registration. Warning! Changes in newMSM: affine method uses SSD by default, Discrete uses Pearson's correlation only. NMI is removed from newMSM, as it was not working in the old version."),
                               false, Utilities::requires_argument);
    Utilities::Option<std::vector<int>> iterations(std::string("--it"), intdefault,
                                        std::string("number of iterations at each resolution (default -–it=3,3,3)"),
                                    false, Utilities::requires_argument);
    std::vector<float> floatdefault;
    Utilities::Option<std::vector<float>> sigma_in(std::string("--sigma_in"),floatdefault,
                                        std::string("smoothing parameter for input image (default --sigma_in=2,2,2)"),
                                    false, Utilities::requires_argument);
    Utilities::Option<std::vector<float>> sigma_ref(std::string("--sigma_ref"),  floatdefault,
                                         std::string("Sigma parameter - smoothing parameter for reference image (set equal to sigma_in by default)"),
                                     false, Utilities::requires_argument);
    Utilities::Option<std::vector<float>> lambda(std::string("--lambda"),  floatdefault,
                                      std::string("Lambda parameter - controls contribution of regulariser "),
                                 false, Utilities::requires_argument);
    Utilities::Option<std::vector<float>> max_dist_pen(std::string("--max_dist_pen"),  floatdefault,
                                                 std::string("Max_distortion_penalty parameter - control weight for the largest distortion "),
                                                 false, Utilities::requires_argument);
    Utilities::Option<std::vector<int>> datagrid(std::string("--datagrid"),intdefault,
                                      std::string("DATA grid resolution (default --datagrid=5,5,5). If parameter = 0 then the native mesh is used."),
                                 false, Utilities::requires_argument);
    Utilities::Option<std::vector<int>> cpgrid(std::string("--CPgrid"),intdefault,
                                    std::string("Control point grid resolution (default --CPgrid=2,3,4)"),
                               false, Utilities::requires_argument);
    Utilities::Option< std::vector<int>> sampgrid(std::string("--SGgrid"),intdefault,
                                  std::string("Sampling grid resolution (default = 2 levels higher than the control point grid)"),
                                  false, Utilities::requires_argument);
    Utilities::Option<std::vector<int>> anatgrid(std::string("--anatgrid"),intdefault,
                                  std::string("Anatomical grid resolution (default = 2 levels higher than the control point grid)"),
                                  false, Utilities::requires_argument);
    std::vector<float> cutthresholddefault(2,0.0); cutthresholddefault[1] = 0.0001;
    Utilities::Option< std::vector<float>> cutthreshold(std::string("--cutthr"),cutthresholddefault,
                                             std::string("Upper and lower thresholds for defining cut vertices (default --cutthr=0,0)"),
                                        false, Utilities::requires_argument);
    Utilities::Option<int> regulariseroption(std::string("--regoption"), 1,
                                  std::string("Choose option for regulariser form lambda*weight*pow(cost,rexp). Where cost can be PAIRWISE or TRI-CLIQUE based. Options are: 1) PAIRWISE - penalising diffences in rotations of neighbouring points (default); 2) TRI_CLIQUE Angle deviation penalty (for spheres); 3) TRI_CLIQUE: Strain-based (for spheres);  4) TRI_CLIQUE Angle deviation penalty (for anatomy); 5) TRI_CLIQUE: Strain-based (for anatomy)"),
                                  false, Utilities::requires_argument);
    Utilities::Option<std::string> doptimizer(std::string("--dopt"),"FastPD",
                                   std::string("discrete optimisation implementation. Choice of: FastPD (default) or HOCR (will reduce to QBPO for pairwise). Warning! ELC and ELC_approx removed from newMSM."),
                              false,Utilities::requires_argument,false);
    Utilities::Option<bool> tricliquelikeihood(std::string("--triclique"), false,
                                    std::string("estimate similarity for triangular patches (rather than circular)"),
                                    false, Utilities::no_argument);
    Utilities::Option<float> shear(std::string("--shearmod"), 0.4,
                        std::string("shear modulus (default 0.4); for use with --regoption 3 "),
                        false,Utilities::requires_argument);
    Utilities::Option<float> bulk(std::string("--bulkmod"), 1.6,
                       std::string("bulk mod (default 1.6); for use with --regoption 3 "),
                       false,Utilities::requires_argument);
    Utilities::Option<float> grouplambda(std::string("--glambda_pairs"), 1,
                              std::string("scaling for pairwise term in group alignment"),
                              false,Utilities::requires_argument,false);
    Utilities::Option<float> kexponent(std::string("--k_exponent"), 2,
                            std::string("exponent inside strain equation (default 2)"),
                            false, Utilities::requires_argument);
    Utilities::Option<float> regulariserexp(std::string("--regexp"), 2.0,
                                 std::string("Regulariser exponent 'rexp' (default 2.0)"),
                                 false,Utilities::requires_argument);
    Utilities::Option<bool> distweight(std::string("--weight"), false,
                            std::string("weight regulariser cost using areal distortion weighting"),
                            false, Utilities::no_argument);
    Utilities::Option<bool> anorm(std::string("--anorm"), false,
                       std::string("norm regulariser cost using mean angle (for HCP compatibility)"),
                       false, Utilities::no_argument);
    Utilities::Option<bool> rescale_labels(std::string("--rescaleL"), false,
                                std::string("rescale label grid rather than using barycentres"),
                                false, Utilities::no_argument);
    Utilities::Option<float> controlptrange(std::string("--cprange"), 1.0,
                                 std::string("Range (as % control point spacing) of data samples (default 1) "),
                                 false,Utilities::requires_argument,false);
    Utilities::Option<bool> intensitynormalize(std::string("--IN"), false,
                                    std::string("Normalize intensity ranges using histogram matching "),
                                    false, Utilities::no_argument);
    Utilities::Option<bool> intensitynormalizewcut(std::string("--INc"), false,
                                        std::string("Normalize intensity ranges using histogram matching excluding cut"),
                                        false, Utilities::no_argument);
    Utilities::Option<bool> variancenormalize(std::string("--VN"), false,
                                   std::string("Variance normalize data "),
                                   false, Utilities::no_argument);
    Utilities::Option<bool> exclude(std::string("--excl"), false,
                         std::string("Ignore the cut when resampling the data"),
                         false, Utilities::no_argument);
    Utilities::Option<bool> quartet(std::string("-Q"), false,
                         std::string("Estimate quartet low rank cost for group reg"),
                         false, Utilities::no_argument,false);
    Utilities::Option<float> affinestepsize(std::string("--stepsize"), 0.01,
                                 std::string("gradient stepping for affine optimisation (default 0.01)"),
                                 false, Utilities::requires_argument);
    Utilities::Option<float> gradsampling(std::string("--gradsampling"), 0.5,
                               std::string("Determines the finite distance spacing for the affine gradient calculation (default 0.5)"),
                               false,Utilities::requires_argument);
    Utilities::Option<std::vector<int>> mciters(std::string("--mciters"), intdefault,
                                   std::string("number of iterations for Monte Carlo optimisation (default == 1000)"),
                                   false,Utilities::requires_argument);
    Utilities::Option<std::vector<float>> labeldist(std::string("--labeldist"), floatdefault,
                                   std::string("distance multiplier for max label distance (default == 0.5)"),
                                   false,Utilities::requires_argument);
    Utilities::Option<int> threads(std::string("--numthreads"), 1,
                        std::string("number of threads for OpenMP (default is single thread)"),
                        false,Utilities::requires_argument);
    //Removed parameters
    Utilities::Option<std::vector<int>> alpha_knn(std::string("--aKNN"),intdefault,
                                    std::string("Warning! This parameter is removed from newMSM."),
                                    false, Utilities::requires_argument);
    Utilities::Option<std::string> meshinterpolationmethod(std::string("--mInt"), "BARY",
                                            std::string("Warning! This parameter is removed from newMSM. (BARY only for surface resampling)"),
                                            false,Utilities::requires_argument);
    Utilities::Option<std::string> datainterpolationmethod(std::string("--dInt"), "ADAP_BARY",
                                            std::string("Warning! This parameter is removed from newMSM. (ADAP_BARY only for data resampling)"),
                                            false,Utilities::requires_argument);
    Utilities::Option<bool> logtransform(std::string("--log"), false,
                              std::string("Warning! This parameter is removed from newMSM. (Not used in code)"),
                               false, Utilities::no_argument);
    Utilities::Option<bool> scaleintensity(std::string("--scale"), false,
                                std::string("Warning! This parameter is removed from newMSM. (Not used in code)"),
                                 false, Utilities::no_argument);
    try {
        // must include all wanted options here (the order determines how
        // the help message is printed)
        options.add(optimizer);
        options.add(simval);
        options.add(iterations);
        options.add(sigma_in);
        options.add(sigma_ref);
        options.add(lambda);
        options.add(max_dist_pen);
        options.add(datagrid);
        options.add(cpgrid);
        options.add(sampgrid);
        options.add(anatgrid);
        options.add(cutthreshold);
        options.add(regulariseroption);
        options.add(doptimizer);
        options.add(tricliquelikeihood);
        options.add(shear);
        options.add(bulk);
        options.add(grouplambda);
        options.add(kexponent);
        options.add(regulariserexp);
        options.add(distweight);
        options.add(anorm);
        options.add(rescale_labels);
        options.add(controlptrange);
        options.add(intensitynormalize);
        options.add(intensitynormalizewcut);
        options.add(variancenormalize);
        options.add(exclude);
        options.add(quartet);
        options.add(affinestepsize);
        options.add(gradsampling);
        options.add(mciters);
        options.add(labeldist);
        options.add(threads);
        //removed parameters
        options.add(alpha_knn);
        options.add(meshinterpolationmethod);
        options.add(datainterpolationmethod);
        options.add(logtransform);
        options.add(scaleintensity);

        if(parameters=="usage")
        {
            options.usage();
            exit(EXIT_SUCCESS); // FSL deploy check
        }

        if (!parameters.empty())
            options.parse_config_file(parameters);
    }
    catch(Utilities::X_OptionError& e)
    {
        options.usage();
        std::cerr << "\n" << e.what() << "\n";
        exit(EXIT_FAILURE);
    }
    catch(std::exception &e)
    {
        std::cerr << "\n" << e.what() << "\n";
        exit(EXIT_FAILURE);
    }

    if(!options.check_compulsory_arguments())
    {
        options.usage();
        exit(EXIT_FAILURE);
    }

    //Removed parameters warning messages.
    if(alpha_knn.set())
        std::cout << "Warning! --aKNN parameter is removed from newMSM." << std::endl;
    if(meshinterpolationmethod.set())
        std::cout << "Warning! --mInt parameter is removed from newMSM. (BARY only for surface resampling)" << std::endl;
    if(datainterpolationmethod.set())
        std::cout << "Warning! --dInt parameter is removed from newMSM. (ADAP_BARY only for surface resampling)" << std::endl;
    if(logtransform.set())
        std::cout << "Warning! --log parameter is removed from newMSM." << std::endl;
    if(scaleintensity.set())
        std::cout << "Warning! --scale parameter is removed from newMSM." << std::endl;

    if(parameters.empty())
    {
        // if no config is supplied use sulc config (as of Sept 2014) as default
        cost.resize(3,"DISCRETE");
        cost.insert(cost.begin(),"RIGID");
        _resolutionlevels=cost.size();
        _lambda.resize(4,0); _lambda[1] = 0.1; _lambda[2] = 0.2; _lambda[3] = 0.3;
        _simval.resize(4,2); _simval[0] = 1;
        _sigma_in.resize(4,2); _sigma_in[2] = 3;_sigma_in[3] = 2;
        _sigma_ref.resize(4,2); _sigma_ref[2] = 1.5; _sigma_ref[3] = 1;
        _iters.resize(4,3); _iters[0] = 50;
        _gridres.resize(4,0); _gridres[1] = 2; _gridres[2] = 3; _gridres[3] = 4;
        _anatres.resize(4,0); _anatres[1] = 4; _anatres[2] = 5; _anatres[3] = 6;
        _genesis.resize(4,4); _genesis[2] = 5; _genesis[3] = 6;
        _sampres.resize(4,0); _sampres[1] = 4; _sampres[2] = 5; _sampres[3] = 6;
    }
    else
    {
        cost = optimizer.value();
        _lambda = lambda.value();
        if(max_dist_pen.set()) max_distortion_penalty = max_dist_pen.value();
        else max_distortion_penalty.resize(_lambda.size(), 0.5);
        // now check for assignments and else set defaults
        if (simval.set()) _simval = simval.value();
        else _simval.resize(cost.size(), 2);
        for(int i = 0; i < _simval.size(); ++i) {
            if (_simval[i] == 3)
                std::cout << "Warning! NMI similarity metric has been removed from newMSM (was not working in old MSM)."
                          << std::endl;
            if (_simval[i] == 1 && cost[i] == "DISCRETE")
                std::cout << "Warning! SSD similarity metric is not supported for discrete method in newMSM (Pearson's correlation only)."
                          << std::endl;
            if (_simval[i] == 2 && (cost[i] == "AFFINE" || cost[i] == "RIGID"))
                std::cout << "Warning! Pearson's correlation similarity metric is not supported for affine/rigid method in newMSM (SSD only)."
                          << std::endl;
        }
        if (iterations.set()) _iters = iterations.value();
        else _iters.resize(cost.size(), 3);
        if (sigma_in.set()) _sigma_in = sigma_in.value();
        else _sigma_in.resize(cost.size(), 2);
        if (sigma_ref.set()) _sigma_ref = sigma_ref.value();
        else _sigma_ref = _sigma_in;
        if (datagrid.set()) _genesis = datagrid.value();
        else _genesis.resize(cost.size(), 5);
        if (cpgrid.set()) _gridres = cpgrid.value();
        else
        {
            _gridres.resize(cost.size(), 2);
            for (unsigned int i = 1; i < cost.size(); i++) _gridres[i] = _gridres[i - 1] + 1;
        }
        if (anatgrid.set()) _anatres = anatgrid.value();
        else
        {
            _anatres.resize(cost.size(), 2);
            for (unsigned int i = 0; i < cost.size(); i++) _anatres[i] = _gridres[i] + 2;
        }
        if (sampgrid.set()) _sampres = sampgrid.value();
        else
        {
            _sampres.resize(cost.size());
            for (unsigned int i = 0; i < cost.size(); i++) _sampres[i] = _gridres[i] + 2;
        }
    }

    _resolutionlevels = cost.size();
    if(grouplambda.set())
    {
        _set_group_lambda=true;
        _pairwiselambda=grouplambda.value();
    }
    _regmode=regulariseroption.value();
    _discreteOPT=doptimizer.value();
    _tricliquelikeihood=tricliquelikeihood.value();
    _shearmod=shear.value();
    _bulkmod=bulk.value();
    _k_exp=kexponent.value();
    if(_discreteOPT=="FastPD") _regmode=1;
    if(intensitynormalizewcut.set())
    {
        _IN=intensitynormalizewcut.value();
        _cut=true;
    }
    else
    {
        _IN=intensitynormalize.value();
        _cut=false;
    }
    _varnorm=variancenormalize.value();
    _exclude=exclude.value();
    _quartet=quartet.value();
    _weight=distweight.value();
    _regoption2norm=anorm.value();
    _threshold=cutthreshold.value();
    _regexp=regulariserexp.value();
    _cprange=controlptrange.value();
    _affinestepsize=affinestepsize.value();
    _affinegradsampling=gradsampling.value();
    _numthreads=threads.value();
    if(mciters.set()) _mciters=mciters.value();
    else _mciters.resize(cost.size(), 50000);
    if(labeldist.set()) _labeldist=labeldist.value();
    else _labeldist.resize(cost.size(), 0.5);
    _rescale_labels=rescale_labels.value();

    if(_verbose)
    {
        std::cout << "\nParameters:\nOptimiser per level: ";
        for(const auto& e : cost) std::cout << e << ' ';
        std::cout << "\nSimilarity measure per level: "; for(const auto& e : _simval) std::cout << e << ' ';
        std::cout << "\nMax iterations per level: "; for(const auto& e : _iters) std::cout << e << ' ';
        std::cout << "\nLambda regulariser: "; for(const auto& e : _lambda) std::cout << e << ' ';
        std::cout << "\nMax distortion penalty: "; for(const auto& e : max_distortion_penalty) std::cout << e << ' ';
        std::cout << "\nSigma for in data per level: "; for(const auto& e : _sigma_in) std::cout << e << ' ';
        std::cout << "\nSigma for ref data per level: "; for(const auto& e : _sigma_ref) std::cout << e << ' ';
        std::cout << "\nDatagrid resolution per level: "; for(const auto& e : _genesis) std::cout << e << ' ';
        std::cout << "\nCP grid resolution per level "; for(const auto& e : _gridres) std::cout << e << ' ';
        std::cout << "\nSampling grid resolution per level: "; for(const auto& e : _sampres) std::cout << e << ' ';
        if (_anat) { std::cout << "\nAnatomical mesh resolution per level: "; for(const auto& e : _anatres) std::cout << e << ' '; }
        if(cutthreshold.set()) { std::cout << "\nCut threshold: "; for(const auto& e : _threshold) std::cout << e << ' '; }
        for(const auto& e : cost) if (e=="AFFINE" || e =="RIGID") {
                std::cout << "\nRigid parameters: Affine step size: " << _affinestepsize
                          << " Affine grad sampling: " << _affinegradsampling;
                break;
            }
        if(variancenormalize.set()) std::cout << "\nVariance normalise set.";
        if(intensitynormalize.set()) std::cout << "\nIntensity normalise set.";
        if(intensitynormalizewcut.set()) std::cout << "\nIntensity normalise with cut set.";
        if(rescale_labels.set()) std::cout << "\nRescale labels set.";
        std::cout << "\nDiscrete implementation: " << _discreteOPT;
        if(_discreteOPT == "MCMC") { std::cout << "\nMonte Carlo iterations: "; for(const auto& e : _mciters) std::cout << e << ' '; }
        std::cout << "\nRegulariser: " <<  _regmode;
        if(_regmode == 3 || _regmode == 5) std::cout << "\nShearmod: " << _shearmod << "; Bulkmod: " << _bulkmod << "; k_exponent: " << _k_exp;
        std::cout << "\nRegulariser exponent: " <<  _regexp;
        std::cout << "\nMultiplier for max label dist: "; for(const auto& e : _labeldist) std::cout << e << ' ';
        std::cout << "\nNumber of execution threads: " << _numthreads << "\n\n";

        std::cout << std::endl;
    }

    if (_regmode > 1 && _discreteOPT == "FastPD")
        throw MeshregException("MeshREG ERROR:: you cannot run higher order clique regularisers with fastPD ");
    if ((int) _threshold.size() != 2)
        throw MeshregException("MeshREG ERROR:: the cut threshold does not contain a limit for upper and lower threshold (too few inputs)");
    if ((int) _simval.size() != _resolutionlevels)
        throw MeshregException("MeshREG ERROR:: config file parameter list lengths are inconsistent: --simval");
    if ((int) _iters.size() != _resolutionlevels)
        throw MeshregException("MeshREG ERROR:: config file  parameter list lengths are inconsistent:--it");
    if ((int) _sigma_in.size() != _resolutionlevels)
        throw MeshregException("MeshREG ERROR:: config file  parameter list lengths are inconsistent: --sigma_in");
    if ((int) _sigma_ref.size() != _resolutionlevels)
        throw MeshregException("MeshREG ERROR:: config file  parameter list lengths are inconsistent:--sigma_ref");
    if ((int) cost.size() != _resolutionlevels)
        throw MeshregException("MeshREG ERROR:: config file parameter list lengths are inconsistent:--opt");
    if ((int) _lambda.size() != _resolutionlevels)
        throw MeshregException("MeshREG ERROR:: config file  parameter list lengths are inconsistent:--lambda");
    if ((int) _genesis.size() != _resolutionlevels)
        throw MeshregException("MeshREG ERROR:: config file  parameter list lengths are inconsistent:--datagrid");
    if ((int) _gridres.size() != _resolutionlevels)
        throw MeshregException("MeshREG ERROR:: config file parameter list lengths are inconsistent:--CPgrid");
    if ((int) _sampres.size() != _resolutionlevels)
        throw MeshregException("MeshREG ERROR:: config file parameter list lengths are inconsistent:--SGres");
}

void Mesh_registration::fix_parameters_for_level(int i) {

    PARAMETERS.clear();

    PARAMETERS.insert(parameterPair("dOPT", _discreteOPT));
    PARAMETERS.insert(parameterPair("lambda", _lambda[i]));
    PARAMETERS.insert(parameterPair("max_penalty", max_distortion_penalty[i]));
    PARAMETERS.insert(parameterPair("lambda_pairs", _pairwiselambda));
    PARAMETERS.insert(parameterPair("set_lambda_pairs", _set_group_lambda));
    PARAMETERS.insert(parameterPair("iters", _iters[i]));
    PARAMETERS.insert(parameterPair("simmeasure", _simval[i]));
    PARAMETERS.insert(parameterPair("sigma_in", _sigma_in[i]));
    PARAMETERS.insert(parameterPair("CPres", _gridres[i]));
    PARAMETERS.insert(parameterPair("SGres", _sampres[i]));
    PARAMETERS.insert(parameterPair("anatres", _anatres[i]));
    PARAMETERS.insert(parameterPair("quartet", _quartet));
    PARAMETERS.insert(parameterPair("regularisermode", _regmode));
    PARAMETERS.insert(parameterPair("TriLikelihood", _tricliquelikeihood));
    PARAMETERS.insert(parameterPair("rescalelabels", _rescale_labels));
    PARAMETERS.insert(parameterPair("shearmodulus", _shearmod));
    PARAMETERS.insert(parameterPair("bulkmodulus", _bulkmod));
    PARAMETERS.insert(parameterPair("range", _cprange));
    PARAMETERS.insert(parameterPair("exponent", _regexp));
    PARAMETERS.insert(parameterPair("weight", _weight));
    PARAMETERS.insert(parameterPair("anorm", _regoption2norm));
    PARAMETERS.insert(parameterPair("scaling", _regscaling));
    PARAMETERS.insert(parameterPair("verbosity", _verbose));
    PARAMETERS.insert(parameterPair("outdir", _outdir));
    PARAMETERS.insert(parameterPair("stepsize", _affinestepsize));
    PARAMETERS.insert(parameterPair("gradsampling", _affinegradsampling));
    PARAMETERS.insert(parameterPair("numthreads", _numthreads));
    PARAMETERS.insert(parameterPair("kexponent", _k_exp));
    PARAMETERS.insert(parameterPair("labeldist", _labeldist[i]));
}

void Mesh_registration::check() {

    if (((MESHES[0].get_coord(0)).norm() - RAD) > 1e-5)
        throw MeshregException("Reg_config ERROR:: input mesh radius has not been normalised to RAD=100");

    if (IN_CFWEIGHTING && (IN_CFWEIGHTING->get_coord(0).norm() - RAD) > 1e-5)
        throw MeshregException("Reg_config ERROR:: input exclusion mesh radius has not been normalised to RAD=100");

    if (REF_CFWEIGHTING && (REF_CFWEIGHTING->get_coord(0).norm() - RAD) > 1e-5)
        throw MeshregException(
                "Reg_config ERROR::reference exclusion mesh radius has not been normalised to RAD=100");
}

void Mesh_registration::set_output_format(const std::string& type) {

    if (type == "GIFTI")
    {
        _surfformat = ".surf.gii";
        _dataformat = ".func.gii";
    }
    else if (type == "ASCII" || type == "ASCII_MAT")
    {
        _surfformat = ".asc";
        if (type == "ASCII")
            _dataformat = ".dpv";
        else
            _dataformat = ".txt"; // for multivariate
    } else
    {
        _surfformat = ".vtk";
        _dataformat = ".txt";
    }
}

// each Matrix can have a different number of rows, for example if the user supplies a
// costfuncion weighting mask for multivariate features only on the reference,
// and also sets exclusion weighting then the reference weighting will be multivariate
// and the source weighting will be univariate - this combines into one mask on the
// source grid (therefore resampling much be reimplemented at every registration step)
NEWMAT::Matrix Mesh_registration::combine_costfunction_weighting(const NEWMAT::Matrix& sourceweight,
                                                          const NEWMAT::Matrix& resampledtargetweight) {

    NEWMAT::Matrix NEW;
    int nrows;
    if (sourceweight.Nrows() >= resampledtargetweight.Nrows())
    {
        NEW = sourceweight;
        nrows = resampledtargetweight.Nrows();
    }
    else
    {
        NEW = resampledtargetweight;
        nrows = sourceweight.Nrows();
    }

    for (int i = 1; i <= NEW.Ncols(); i++)
        for (int j = 1; j <= nrows; j++)
            NEW(j, i) = (sourceweight(j, i) + resampledtargetweight(j, i)) / 2.0;

    return NEW;
}

std::vector<std::string> Mesh_registration::read_ascii_list(const std::string& filename) {

    std::vector<std::string> list;

    std::ifstream fs(filename);
    std::string tmp;

    if(fs)
    {
        fs >> tmp;
        do
        {
            list.push_back(tmp);
            fs >> tmp;
        }
        while(!fs.eof());
    }
    else
        throw MeshregException("Error reading ascii file");

    return list;
}

} //namespace newmeshreg
