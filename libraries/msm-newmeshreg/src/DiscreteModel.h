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
#ifndef NEWMESHREG_DISCRETEMODEL_H
#define NEWMESHREG_DISCRETEMODEL_H

#include <omp.h>

#include "DiscreteCostFunction.h"

namespace newmeshreg {

class DiscreteModel {

public:
    virtual ~DiscreteModel() {
        delete[] labeling;
        delete[] pairs;
        delete[] triplets;
    }

    //---GET---//
    int getNumNodes() const { return m_num_nodes; }
    int getNumLabels() const { return m_num_labels; }
    int getNumPairs() const { return m_num_pairs; }
    int getNumTriplets() const { return m_num_triplets; }

    int* getLabeling() { return labeling; }
    const int* getLabeling() const { return labeling; }
    const int* getPairs() const { return pairs; }
    const int* getTriplets() const { return triplets; }

    virtual inline std::shared_ptr<DiscreteCostFunction> getCostFunction() = 0;

    //---COMPUTE---//
    virtual inline void computeUnaryCosts() {}
    virtual inline double computeUnaryCost(int node, int label) { return 0; }
    virtual inline void computePairwiseCosts() {}
    virtual inline double computePairwiseCost(int pair, int labelA, int labelB) { return 0; }
    virtual inline void computeTripletCosts() {}
    virtual inline double computeTripletCost(int triplet, int labelA, int labelB, int labelC) { return 0; }
    virtual double inline computeTripletCostTri(int trID, int labelA, int labelB, int labelC) { return 0; }
    virtual inline double evaluateTotalCostSum() { return 0; }

    //---MODIFY---//
    virtual void applyLabeling(int *discreteLabeling) {}

protected:
    void resetLabeling() { if(labeling) std::fill(labeling,labeling + m_num_nodes,0.0); }
    void initLabeling() {
        if(m_num_nodes != 0)
        {
            delete[] labeling;
            labeling = new int[m_num_nodes];
            resetLabeling();
        }
    }

    int m_num_nodes = 0;    // Number of model nodes (e.g. grid nodes, variables, etc).
    int m_num_labels = 0;   // Number of labels.
    int m_num_pairs = 0;    // Number of node pairs.
    int m_num_triplets = 0; // Number of node triplets.
    int* labeling = nullptr;      // Labeling array.
    int* pairs = nullptr;         // Node pairs array.
    int* triplets = nullptr;      // Node triplets array.

    std::string m_outdir;
    bool m_verbosity = false;
    int _nthreads = 1;
};

class NonLinearSRegDiscreteModel : public DiscreteModel {

public:
    NonLinearSRegDiscreteModel() = default;
    explicit NonLinearSRegDiscreteModel(myparam& PAR) {
        set_parameters(PAR);
        initialize_cost_function(m_multivariate, PAR);
    }

    virtual void Initialize(const newresampler::Mesh&);
    void initialize_cost_function(bool MV, myparam &P);
    void set_parameters(myparam& PAR);
    void Initialize_sampling_grid();
    void label_sampling_grid(int, double, newresampler::Mesh&);
    std::vector<newresampler::Point> rescale_sampling_grid();

    virtual void setupCostFunction();

    //---COMPUTE COSTS---//
    void inline computeUnaryCosts() override { costfct->computeUnaryCosts(); }
    double inline computeUnaryCost(int node, int label) override { return costfct->computeUnaryCost(node, label); }
    void inline computePairwiseCosts() override { costfct->computePairwiseCosts(pairs); }
    double inline computePairwiseCost(int pair, int labelA, int labelB) override { return costfct->computePairwiseCost(pair, labelA, labelB); }
    void inline computeTripletCosts() override { costfct->computeTripletCosts(); }
    double inline computeTripletCost(int triplet, int labelA, int labelB, int labelC) override { return costfct->computeTripletCost(triplet, labelA, labelB, labelC); }
    double inline evaluateTotalCostSum() override { return costfct->evaluateTotalCostSum(labeling, pairs, triplets); }

    virtual inline void set_meshspace(const newresampler::Mesh& target, const newresampler::Mesh& source, int num) { m_TARGET = target; m_SOURCE = source; }
    inline void set_anatomical_meshspace(const newresampler::Mesh& ref_sphere, const newresampler::Mesh& ref_anat,
                                  const newresampler::Mesh& source_sphere, const newresampler::Mesh & source_anat) {
        costfct->set_anatomical(ref_sphere, ref_anat, source_sphere, source_anat);
    }
    inline void set_anatomical_neighbourhood(const std::vector<std::map<int,double>>& weights, const std::vector<std::vector<int>>& neighbourhood) {
        costfct->set_anatomical_neighbourhood(weights, neighbourhood);
    }
    inline void set_featurespace(const std::shared_ptr<featurespace>& FEATURES) {
        FEAT = FEATURES;
        costfct->set_featurespace(FEATURES);
    }
    // costfunction weighting combines source and reference weightings at beginning of optimisation iteration -
    // will not be 100% accurate but will remove any sensitivity of the label choices to weighting
    inline void setupCostFunctionWeighting(const NEWMAT::Matrix& Weight) { costfct->set_dataaffintyweighting(Weight); }

    // source needs to be reset after every iteration of discrete optimisation
    virtual inline void reset_meshspace(const newresampler::Mesh& source, int num) {
        m_SOURCE = source;
        costfct->reset_source(source,num);
    }

    virtual inline void reset_CPgrid(const newresampler::Mesh& grid, int num) { m_CPgrid = grid; }
    virtual inline void warp_CPgrid(newresampler::Mesh& START, newresampler::Mesh& END, int num) {
        newresampler::sphere_project_warp(m_CPgrid,START,END, _nthreads);
        unfold(m_CPgrid, m_verbosity);
    }
    inline void set_debug(){ m_debug = true; } // for debuging
    inline newresampler::Mesh get_TARGET() { return m_TARGET; }
    virtual inline newresampler::Mesh get_CPgrid(int num) { return m_CPgrid; }
    virtual void set_masks(const newresampler::Mesh& m) {}
    inline std::shared_ptr<DiscreteCostFunction> getCostFunction() override { return costfct; }
    virtual void applyLabeling();

    double getMCParam() const { return MCparam; }

protected:
    newresampler::Mesh m_TARGET; // TARGET MESH
    newresampler::Mesh m_SOURCE; // SOURCE MESH
    newresampler::Mesh m_CPgrid; // CONTROL POINT GRID
    newresampler::Mesh m_samplinggrid;

    std::shared_ptr<newresampler::Octree> m_inputtree;
    std::shared_ptr<featurespace> FEAT;

    int m_SGres = 4; // sampling grid resolution
    int m_iter = 0; // iteration of the discrete optimisation
    int m_centroid = 0; // used for selecting which sampling grid vertex will form the center of the sampling grid
    int m_regoption = 2; // sim measure i.e. correlation
    double m_maxs_dist = 0.0; // define maximum distance between the centre of the sampling grid and the furthest label
    double _labeldist = 0.5;
    double range = 1;
    double m_scale = 0.0;
    double MCparam = 0.8;
    bool m_multivariate = false;
    bool m_debug = false;
    bool m_triclique = false;
    std::string optimiser;
    bool _pairwise = false;
    bool m_rescalelabels = false;
    newresampler::Point centre;

    std::vector<newresampler::Point> m_samples;  // samples based on  vertices of sampling grid
    std::vector<newresampler::Point> m_barycentres; // samples based on barycentres of sampling grid
    std::vector<newresampler::Point> m_labels; // labels iterates between samples and barycnetres and is the label set used within cosfct
    std::vector<NEWMAT::Matrix> m_ROT; // rotates sampling grid to each control point
    std::shared_ptr<NonLinearSRegDiscreteCostFunction> costfct;  // costfunction object

    //---INIT---//
    virtual void estimate_pairs();
    virtual void estimate_triplets();
    virtual void get_rotations();
};

} //namespace newmeshreg

#endif //NEWMESHREG_DISCRETEMODEL_H
