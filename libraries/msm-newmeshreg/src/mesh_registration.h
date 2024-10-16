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
#ifndef NEWMESHREG_MESHREG_H
#define NEWMESHREG_MESHREG_H

#include <utils/options.h>

#include "rigid_costfunction.h"
#include "newresampler/resampler.h"
#include "DiscreteModel.h"
#include "mcmc_opt.h"

#ifdef HAS_HOCR
#include "Fusion/Fusion.h"
#endif
#ifdef HAS_FPD
#include "FastPD/FastPD.h"
#endif

namespace newmeshreg {

class Mesh_registration {

public:
    Mesh_registration();

    //---ENTRY POINT---//
    void run_multiresolutions(const std::string& parameters);

    //---SET FUNCTIONS---//
    void set_input(const newresampler::Mesh& M);
    void set_input(const std::string& M);
    void set_reference(const newresampler::Mesh& M);
    void set_reference(const std::string& M);
    void set_anatomical(const std::string& M1, const std::string& M2);
    void set_transformed(const std::string& M);
    void set_input_cfweighting(const std::string& E);
    void set_reference_cfweighting(const std::string& E);
    void set_output_format(const std::string& type);
    inline void print_config_options() { parse_reg_options("usage"); }
    inline void set_verbosity(bool V) { _verbose = V; }
    inline void is_sparse(bool sp) { _issparse = sp; } // input data is sparse
    inline void set_outdir(const std::string& s) { _outdir = s; }
    inline void set_debug(bool test) { _debug = test; }
    inline void set_CMpathin(const std::string& s) { CMfile_in = s; }
    inline void set_CMpathref(const std::string& s) { CMfile_ref = s; }

protected:
    int level = 0;
    bool isrigid = false;
    std::shared_ptr<NonLinearSRegDiscreteModel> model;
    std::shared_ptr<Rigid_cost_function> rigidcf;

    //---DATA AND MESHES---//
    std::vector<newresampler::Mesh> MESHES;  // original input (moving) mesh

    newresampler::Mesh transformed_mesh;  // original input (moving) mesh in transformed position (i.e. from previous alignment step)
    newresampler::Mesh in_anat; // input anatomical mesh
    newresampler::Mesh ref_anat; // reference anatomical mesh
    newresampler::Mesh SPH_orig; // original low res icosphere mesh
    newresampler::Mesh SPH_reg;  // transformed low res icosphere mesh
    newresampler::Mesh ANAT_orig; // original low res icosphere mesh

    std::shared_ptr<newresampler::Mesh> IN_CFWEIGHTING; // cost function weights for high resolution meshes
    std::shared_ptr<newresampler::Mesh> REF_CFWEIGHTING;
    NEWMAT::Matrix SPHin_CFWEIGHTING; // downsampled costfunction weightings
    NEWMAT::Matrix SPHref_CFWEIGHTING;

    std::shared_ptr<featurespace> FEAT;
    std::vector<std::string> DATAlist;

    std::string CMfile_in; // location of connectivity/data matrix for input
    std::string CMfile_ref; // location of connectivity/data matrix for reference
    std::string _outdir;

    std::string _surfformat = ".surf";
    std::string _dataformat = ".func";

    //---DATA OPTIONS---//
    bool _IN = false;
    bool _issparse = false;

    std::string _discreteOPT;

    //---REGISTRATION PARAMETERS---//
    myparam PARAMETERS;
    std::vector<std::string> cost;   // controls registration method i.e. rigid, discrete
    std::vector<int> _genesis;           // ico mesh resolution at this level
    std::vector<float> _sigma_in;        // smoothing of input
    std::vector<float> _sigma_ref;      // smoothing of reference
    std::vector<int> _simval; // code determines how similarity is assessed 1=SSD; 2=Pearson's correlation;
    std::vector<float> _lambda;         // controls regularisation
    std::vector<float> _threshold;         // controls cut exclusion (2D upper and lower thresholds for defining cut vertices)
    std::vector<int> _iters; // total per resolution level
    std::vector<int> _gridres; // control point grid resolution (for discrete reg)
    std::vector<int> _anatres; // control point grid resolution (for discrete reg)
    std::vector<int> _sampres; // sampling grid for discrete reg (should be higher than control point)
    bool _exclude = false;  // exclusion zone
    bool _cut = false;  // exclusion zone
    bool _varnorm = false; // variance normalise
    bool _tricliquelikeihood = false;
    bool _patchwise = false;
    bool _anat = false; // if anatomical meshes provided
    bool _rescale_labels = false;
    bool _incfw = false;    //if cost function weightings provided
    bool _refcfw = false;
    bool fixnan = false;
    double _k_exp = 2.0;
    double percentile = 1.0;
    int _resolutionlevels = 1;  // default 1
    int _regmode = 0; //regulariser option

    //---REGULARISER OPTIONS---//
    double _regexp = 0.0; // choice of exponent
    double _regscaling = 0.0; // choice of exponent scaling for options 2 and 3
    double _shearmod = 0.4; // for strain regulariser
    double _bulkmod = 1.6;    // for strain regulariser
    double _cprange = 1.0;
    std::vector<int> _mciters;  // Monte Carlo iterations
    double _mcparameter = 0.8;  // param for Monte Carlo random distribution

    //---AFFINE PARAMETERS---//
    double _affinestepsize = 0.0;
    double _affinegradsampling = 0.0;

    //---MISC---//
    int _numthreads = 1;
    bool _verbose = false;
    bool _debug = false;

    //---INITIALISATION---//
    void parse_reg_options(const std::string& parameters);
    void fix_parameters_for_level(int i);
    void check(); // checks you have all the necessary data
    std::vector<std::string> read_ascii_list(const std::string& filename);

    //---PREPROCESSING---//
    virtual void initialize_level(int current_lvl);
    NEWMAT::Matrix combine_costfunction_weighting(const NEWMAT::Matrix &, const NEWMAT::Matrix &);
    newresampler::Mesh resample_anatomy(const newresampler::Mesh& control_grid, std::vector<std::map<int,double>>& baryweights, std::vector<std::vector<int>>& ANAT_to_CPgrid_neighbours, int current_lvl);
    NEWMAT::Matrix downsample_cfweighting(const newresampler::Mesh& SPH, std::shared_ptr<newresampler::Mesh> CFWEIGHTING, std::shared_ptr<newresampler::Mesh> EXCL);
    newresampler::Mesh project_CPgrid(newresampler::Mesh SPH_in, const newresampler::Mesh& REG, int num = 0);
    NEWMAT::Matrix combine_weighting();

    //---RUN---//
    virtual void evaluate();
    virtual void transform(const std::string& filename);
    virtual void run_discrete_opt();

    //---POSTPROCESSING---//
    virtual inline void saveSPH_reg(const std::string& filename) const { SPH_reg.save(filename + "sphere.LR.reg" + _surfformat); }
    virtual void save_transformed_data(const std::string& filename);
};

} //namespace newmeshreg

#endif //NEWMESHREG_MESHREG_H
