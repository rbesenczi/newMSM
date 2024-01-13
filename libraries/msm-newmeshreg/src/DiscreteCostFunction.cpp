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
#include "DiscreteCostFunction.h"

namespace newmeshreg {

void DiscreteCostFunction::initialize(int numNodes, int numLabels, int numPairs, int numTriplets) {

    if (m_num_nodes != numNodes || m_num_labels != numLabels)
    {
        delete[] unarycosts;
        unarycosts = new double[numNodes * numLabels];
    }

    if (m_num_pairs != numPairs || m_num_labels != numLabels)
    {
        delete[] paircosts;
        paircosts = new double[numPairs * numLabels * numLabels];
    }

    if (m_num_triplets != numTriplets || m_num_labels != numLabels)
    {
        delete[] tripletcosts;
        tripletcosts = new double[numTriplets * numLabels * numLabels * numLabels];
    }

    m_num_nodes = numNodes;
    m_num_labels = numLabels;
    m_num_pairs = numPairs;
    m_num_triplets = numTriplets;

    std::fill(unarycosts,unarycosts+m_num_labels*m_num_nodes,0.0f);
    std::fill(paircosts,paircosts+m_num_labels*m_num_labels*m_num_pairs,0.0f);
    std::fill(tripletcosts,tripletcosts+m_num_triplets*m_num_labels*m_num_labels*m_num_labels,0.0f);
}

void DiscreteCostFunction::reset() {
    if(unarycosts) std::fill(unarycosts,unarycosts+m_num_labels*m_num_nodes,0.0f);
    if(paircosts) std::fill(paircosts,paircosts+m_num_labels*m_num_labels*m_num_pairs,0.0f);
    if(tripletcosts) std::fill(tripletcosts,tripletcosts+m_num_triplets*m_num_labels*m_num_labels*m_num_labels,0.0f);
}

double DiscreteCostFunction::evaluateTotalCostSum(const int *labeling, const int *pairs, const int *triplets) {

    double cost_sum_unary = 0.0f;
    double cost_sum_pairwise = 0.0f;
    double cost_sum_triplet = 0.0f;

    for (int i = 0; i < m_num_nodes; ++i) {
        double unary_cost = computeUnaryCost(i, labeling[i]);
        cost_sum_unary += unary_cost;
    }

    for (int p = 0; p < m_num_pairs; ++p) {
        double pair_cost = computePairwiseCost(p, labeling[pairs[p * 2]], labeling[pairs[p * 2 + 1]]);
        cost_sum_pairwise += pair_cost;
    }

    for (int t = 0; t < m_num_triplets; ++t) {
        double triplet_cost = computeTripletCost(t, labeling[triplets[t * 3]], labeling[triplets[t * 3 + 1]],
                                                 labeling[triplets[t * 3 + 2]]);
        cost_sum_triplet += triplet_cost;
    }

    if(_verbosity)
        std::cout << "cost_sum_unary " << cost_sum_unary << " cost_sum_pairwise " << cost_sum_pairwise
        << " cost_sum_triplet " << cost_sum_triplet
        <<" total " <<  cost_sum_unary + cost_sum_pairwise + cost_sum_triplet
        << " m_num_triplets 2 " << m_num_triplets <<  std::endl;

    return cost_sum_unary + cost_sum_pairwise + cost_sum_triplet;
}

newresampler::Mesh NonLinearSRegDiscreteCostFunction::project_anatomical() {
    newresampler::Mesh _aICOtrans = _aICO;
    barycentric_mesh_interpolation(_aICOtrans,_ORIG,_SOURCE, _threads);
    return newresampler::project_mesh(_aICOtrans,_TARGEThi,_aTARGET, _threads);
}

void NonLinearSRegDiscreteCostFunction::reset_anatomical() {

    if(_aSOURCE.nvertices() > 0)
    {
        NEWMAT::ColumnVector strainstmp(_aSOURCE.ntriangles());
        _aSOURCEtrans = project_anatomical();
        MAXstrain = 0.0;

        for (int i = 0; i < _aSOURCE.ntriangles(); i++)
            strainstmp(i + 1) = calculate_triangular_strain(i, _aSOURCE, _aSOURCEtrans, _mu, _kappa);

        for (int i = 0; i < _aSOURCE.ntriangles(); i++)
            if(strainstmp(i+1) > MAXstrain)
                MAXstrain = strainstmp(i+1);
    }
}

void NonLinearSRegDiscreteCostFunction::set_initial_angles(const std::vector<std::vector<double>>& angles) {
    double meanang = 0.0;
    for (const auto& angle : angles)
        for (double j : angle)
            meanang += j;

    _MEANANGLE = meanang / (3.0*angles.size());
}

bool NonLinearSRegDiscreteCostFunction::within_controlpt_range(int CPindex, int sourceindex) {

    return ((2 * RAD *
             asin((_CPgrid.get_coord(CPindex)-_SOURCE.get_coord(sourceindex)).norm() / (2*RAD)))
             < _controlptrange * MAXSEP(CPindex + 1));
}

void NonLinearSRegDiscreteCostFunction::initialize(int numNodes, int numLabels, int numPairs, int numTriplets) {
    if (_TARGET.nvertices() == 0 || _SOURCE.nvertices() == 0)
        throw MeshregException("CostFunction::You must supply source and target meshes.");
    if(_HIGHREScfweight.Ncols() != _SOURCE.nvertices())
        _HIGHREScfweight.ReSize(1, _SOURCE.nvertices());
    if (_HIGHREScfweight.Nrows() != 1 && _HIGHREScfweight.Nrows() != FEAT->get_dim())
        throw MeshregException("DiscreteModel ERROR:: costfunction weighting has dimensions incompatible with data");
    DiscreteCostFunction::initialize(numNodes, numLabels, numPairs, numTriplets);
}

void NonLinearSRegDiscreteCostFunction::set_parameters(myparam& ALLPARAMS) {
    myparam::iterator it;
    it=ALLPARAMS.find("exponent");_rexp=std::get<float>(it->second);
    it=ALLPARAMS.find("weight");_dweight=std::get<bool>(it->second);
    it=ALLPARAMS.find("anorm");_anorm=std::get<bool>(it->second);
    it=ALLPARAMS.find("shearmodulus");_mu=std::get<float>(it->second);
    it=ALLPARAMS.find("bulkmodulus");_kappa=std::get<float>(it->second);
    it=ALLPARAMS.find("kexponent");_k_exp=std::get<float>(it->second);
    it=ALLPARAMS.find("lambda"); _reglambda=std::get<float>(it->second);
    it=ALLPARAMS.find("range"); _controlptrange=std::get<float>(it->second);
    it=ALLPARAMS.find("simmeasure"); _simmeasure=std::get<int>(it->second); sim.set_simval(_simmeasure);
    it=ALLPARAMS.find("verbosity"); _verbosity=std::get<bool>(it->second);
    it=ALLPARAMS.find("regularisermode"); _rmode=std::get<int>(it->second);
    it=ALLPARAMS.find("numthreads"); _threads=std::get<int>(it->second);
    it=ALLPARAMS.find("dOPT"); dopt=std::get<std::string>(it->second);
}

double NonLinearSRegDiscreteCostFunction::computeTripletCost(int triplet, int labelA, int labelB, int labelC) {

    double cost = 0.0, weight = 1.0;
    if(triplet == 0 && _debug) { sumlikelihood = 0.0; sumregcost = 0.0; }

    std::map<int,newresampler::Point> vertex;

    std::vector<int> id = { _triplets[3*triplet], _triplets[3*triplet+1], _triplets[3*triplet+2] };

    newresampler::Point v0 = _CPgrid.get_coord(id[0]);
    newresampler::Point v1 = _CPgrid.get_coord(id[1]);
    newresampler::Point v2 = _CPgrid.get_coord(id[2]);

    vertex[id[0]] = (*ROTATIONS)[id[0]] * _labels[labelA];
    vertex[id[1]] = (*ROTATIONS)[id[1]] * _labels[labelB];
    vertex[id[2]] = (*ROTATIONS)[id[2]] * _labels[labelC];

    newresampler::Triangle TRI(vertex[id[0]],
                               vertex[id[1]],
                               vertex[id[2]], 0);
    newresampler::Triangle TRI_noDEF(v0,v1,v2,0);

    double likelihood = triplet_likelihood(triplet, id[0], id[1], id[2], vertex[id[0]], vertex[id[1]], vertex[id[2]]);

    // only estimate cost if it doesn't cause folding
    if((TRI.normal() | TRI_noDEF.normal()) < 0)
    {
        cost = 1.0;
        weight += 1e6;
    }
    else
    {
        switch(_rmode)
        {
            case 2:
            {
                double tweight = 1.0, diff = 0.0;
                newresampler::Triangle TRI_ORIG(_ORIG.get_coord(_triplets[3*triplet]),
                                                _ORIG.get_coord(_triplets[3*triplet+1]),
                                                _ORIG.get_coord(_triplets[3*triplet+2]), 0);
                std::vector<double> deformed_angles = TRI.get_angles();
                std::vector<double> orig_angles = TRI_ORIG.get_angles();

                double distortion = log2(TRI.get_area()/TRI_ORIG.get_area());
                if (_dweight && distortion != 1) tweight = exp(abs(distortion - 1));
                else tweight = 1.0;

                for (int i = 0; i < 3; i++)
                    diff += std::pow(deformed_angles[i] - orig_angles[i], 2);

                cost = tweight*(sqrt(diff));

                if (_anorm) weight *= (1 / _MEANANGLE);
                break;
            }
            case 3:
            {
                newresampler::Triangle TRI_ORIG(_ORIG.get_coord(_triplets[3*triplet]),
                                                _ORIG.get_coord(_triplets[3*triplet+1]),
                                                _ORIG.get_coord(_triplets[3*triplet+2]), 0);

                cost = calculate_triangular_strain(TRI_ORIG, TRI, _mu, _kappa, std::shared_ptr<NEWMAT::ColumnVector>(), _k_exp);
                break;
            }
            case 4:
            {
                double tweight = 1.0, diff = 0.0;
                std::map<int,bool> moved2;
                std::map<int,newresampler::Point> transformed_points;

                for (unsigned int n = 0; n < NEARESTFACES[triplet].size(); n++)
                {
                    newresampler::Triangle TRIorig = _aSOURCE.get_triangle(NEARESTFACES[triplet][n]);
                    newresampler::Triangle TRItrans = deform_anatomy(triplet, n, vertex, moved2, transformed_points);
                    std::vector<double> orig_angles = TRIorig.get_angles();
                    std::vector<double> deformed_angles = TRItrans.get_angles();

                    double distortion = log2(TRItrans.get_area()/TRIorig.get_area());
                    if (_dweight && distortion != 1) tweight = exp(abs(distortion - 1));
                    else tweight = 1;

                    for (int i = 0; i < 3; i++)
                        diff += std::pow(deformed_angles[i] - orig_angles[i], 2);
                    diff = tweight * sqrt(diff);
                    cost += diff;   // TODO this is a bit odd, as diff is not zeroed at any point in this loop. Maybe double summing the cost?
                }
                cost = cost / (double) NEARESTFACES[triplet].size();
                break;
            }
            case 5:
            {
                std::map<int,bool> moved2;
                std::map<int,newresampler::Point> transformed_points;

                for (unsigned int n = 0; n < NEARESTFACES[triplet].size(); n++)
                {
                    newresampler::Triangle TRIorig = _aSOURCE.get_triangle(NEARESTFACES[triplet][n]);
                    newresampler::Triangle TRItrans = deform_anatomy(triplet, n, vertex, moved2, transformed_points);
                    cost += calculate_triangular_strain(TRIorig, TRItrans, _mu, _kappa, std::shared_ptr<NEWMAT::ColumnVector>(), _k_exp);
                }
                cost = cost / (double) NEARESTFACES[triplet].size();
                break;
            }
            default:
                throw MeshregException("DiscreteModel computeTripletCost regoption does not exist");
        }
    }

    if (abs(cost) < 1e-8) cost = 0.0;
    if(_debug)
    {
        sumlikelihood += likelihood;
        sumregcost += weight * _reglambda * MISCMATHS::pow(cost,_rexp);
    }

    return likelihood + weight * _reglambda * MISCMATHS::pow(cost,_rexp); // normalise to try and ensure equivalent lambda for each resolution level
}

double NonLinearSRegDiscreteCostFunction::computePairwiseCost(int pair, int labelA, int labelB) {

    NEWMAT::Matrix R1 = estimate_rotation_matrix(_CPgrid.get_coord(_pairs[2*pair]),(*ROTATIONS)[_pairs[2*pair]]*_labels[labelA]);
    NEWMAT::Matrix R2 = estimate_rotation_matrix(_CPgrid.get_coord(_pairs[2*pair+1]),(*ROTATIONS)[_pairs[2*pair+1]]*_labels[labelB]);
    NEWMAT::Matrix R_diff = R1.t() * R2;

    const double theta_MVD = 2 * asin(MVDmax/(2*RAD));
    const double theta = acos((R_diff.Trace()-1)/2);
    double weight = 1.0, totalcost = 0.0, cost = 0.0;

    newresampler::Point v0 = _CPgrid.get_coord(_pairs[2*pair]);
    newresampler::Point v1 = _CPgrid.get_coord(_pairs[2*pair+1]);

    _CPgrid.set_coord(_pairs[2*pair],(*ROTATIONS)[_pairs[2*pair]]*_labels[labelA]);
    _CPgrid.set_coord(_pairs[2*pair+1],(*ROTATIONS)[_pairs[2*pair+1]]*_labels[labelB]);

    if(fabs(1 - (R_diff.Trace() - 1) / 2) > EPSILON)
    {
        for(auto j = _CPgrid.tIDbegin(_pairs[2*pair]); j != _CPgrid.tIDend(_pairs[2*pair]); j++)
        {
            if((_oCPgrid.get_triangle(*j).normal() | _CPgrid.get_triangle(*j).normal()) < 0)
            {
                cost = 1.0;
                weight += 1e6;
                break;
            }
        }

        cost = (sqrt(2) * theta) / theta_MVD;

        if(_rexp == 1)
            totalcost = weight * _reglambda * cost;
        else
            totalcost = weight * _reglambda * MISCMATHS::pow(cost,_rexp);
    }

    _CPgrid.set_coord(_pairs[2*pair],v0);
    _CPgrid.set_coord(_pairs[2*pair+1],v1);

    return totalcost;
}

void NonLinearSRegDiscreteCostFunction::computePairwiseCosts(const int *pairs) {

    for (int i = 0; i < m_num_pairs; i++)
        for (unsigned int j = 0; j < _labels.size(); j++)
            for (unsigned int k = 0; k < _labels.size(); k++)
                paircosts[i * m_num_labels * m_num_labels + k * m_num_labels + j] = computePairwiseCost(i, j, k);
}

void NonLinearSRegDiscreteCostFunction::computeUnaryCosts() {
    // for each control point resample data into blocks each influencing a single control point
    // calculates similarity
    for (int j = 0; j < m_num_labels; j++)
        #pragma omp parallel for num_threads(_threads)
        for (int k = 0; k < _CPgrid.nvertices(); k++)
            unarycosts[j * m_num_nodes + k] = computeUnaryCost(k, j);
}

newresampler::Triangle NonLinearSRegDiscreteCostFunction::deform_anatomy(int trip, int n, std::map<int,newresampler::Point>& vertex,
                                                                         std::map<int,bool>& moved, std::map<int,newresampler::Point>& transformed) {
    newresampler::Triangle TRItrans;

    for(int i=0;i<3;i++)
    { // for each point in face
        const int tindex = _aSOURCE.get_triangle(NEARESTFACES[trip][n]).get_vertex_no(i);

        if(moved.find(tindex) == moved.end())
        {
            moved[tindex] = true;
            newresampler::Point newPt;

            //newPt.X = 0; newPt.Y = 0; newPt.Z = 0;
            for (const auto& it: _ANATbaryweights[tindex])
                newPt += vertex[it.first] * it.second;

            newresampler::Triangle closest_triangle;
            try {
                closest_triangle = anattree->get_closest_triangle(newPt);
            } catch (newresampler::MeshException& e) {
                std::cout << "Warning! Cannot find closest triangle on anatomical mesh. This is a known bug in the Octree search algorithm in anatomical MSM."
                      << "\n\tUsing default Triangle with ID==0." << std::endl;
                closest_triangle.set(newresampler::Point{0,0,0},newresampler::Point{0,0,0},newresampler::Point{0,0,0},0);
            }

            newresampler::Point v0 = closest_triangle.get_vertex_coord(0),
                                v1 = closest_triangle.get_vertex_coord(1),
                                v2 = closest_triangle.get_vertex_coord(2);
            int n0 = closest_triangle.get_vertex_no(0),
                n1 = closest_triangle.get_vertex_no(1),
                n2 = closest_triangle.get_vertex_no(2);

            std::map<int,double> weight = newresampler::calc_barycentric_weights(v0,v1,v2,newPt,n0,n1,n2);

            newPt.X = 0; newPt.Y = 0; newPt.Z = 0;
            for (const auto& it: weight)
                newPt += _aTARGET.get_coord(it.first) * it.second;

            transformed[tindex] = newPt;
            TRItrans.set_vertex(i, newPt);
        }
        else
            TRItrans.set_vertex(i, transformed[tindex]);
    }

    return TRItrans;
}

void NonLinearSRegDiscreteCostFunction::resample_weights(){

    // TAKE DATA TERM WEIGHTING AS MAXIMUM OF WEIGHTS ACROSS ALL DIMENSIONS
    AbsoluteWeights.ReSize(_SOURCE.nvertices());
    AbsoluteWeights = 0;

    for (int k = 1; k <= _SOURCE.nvertices(); k++)
    {
        double maxweight = std::numeric_limits<double>::lowest();
        for (int j = 1; j <= _HIGHREScfweight.Nrows(); j++)
            if (_HIGHREScfweight(j, k) > maxweight)
                maxweight = _HIGHREScfweight(j, k);

        AbsoluteWeights(k) = maxweight;
    }

    newresampler::Mesh tmp = _SOURCE;
    tmp.set_pvalues(AbsoluteWeights);
    AbsoluteWeights = newresampler::metric_resample(tmp, _CPgrid, _threads).get_pvalues();
}

void UnivariateNonLinearSRegDiscreteCostFunction::initialize(int numNodes, int numLabels, int numPairs, int numTriplets) {

    NonLinearSRegDiscreteCostFunction::initialize(numNodes, numLabels, numPairs, numTriplets);
    _sourcedata.clear(); _sourcedata.resize(_CPgrid.nvertices());
    _sourceinrange.clear(); _sourceinrange.resize(_CPgrid.nvertices());
    _targetdata.clear(); _targetdata.resize(_CPgrid.nvertices());
    _weights.clear(); _weights.resize(_CPgrid.nvertices());
}

void UnivariateNonLinearSRegDiscreteCostFunction::get_source_data() {

    for (auto& i : _sourcedata) i.clear();

    for (int k = 0; k < _CPgrid.nvertices(); k++)
        for (int i = 0; i < _SOURCE.nvertices(); i++)
            if (within_controlpt_range(k, i))
            {
                _sourceinrange[k].push_back(i);
                _sourcedata[k].push_back(FEAT->get_input_val(1, i + 1));
                if (_HIGHREScfweight.Nrows() >= 1)
                    _weights[k].emplace_back(_HIGHREScfweight(1, i + 1));
                else
                    _weights[k].emplace_back(1.0);
            }
    resample_weights();
}

void UnivariateNonLinearSRegDiscreteCostFunction::get_target_data(int node, const NEWMAT::Matrix& PtROTATOR) {

    _targetdata[node].clear();
    _targetdata[node].resize(_sourceinrange[node].size());

    for(unsigned int i = 0; i < _sourceinrange[node].size(); i++)
    {
        newresampler::Point tmp = PtROTATOR * _SOURCE.get_coord(_sourceinrange[node][i]);

        newresampler::Triangle closest_triangle = targettree->get_closest_triangle(tmp);

        newresampler::Point v0 = closest_triangle.get_vertex_coord(0),
                            v1 = closest_triangle.get_vertex_coord(1),
                            v2 = closest_triangle.get_vertex_coord(2);
                        int n0 = closest_triangle.get_vertex_no(0),
                            n1 = closest_triangle.get_vertex_no(1),
                            n2 = closest_triangle.get_vertex_no(2);

        _targetdata[node][i] = newresampler::barycentric_interpolation(v0, v1, v2, tmp,
                                                                FEAT->get_ref_val(1, n0+1),
                                                                FEAT->get_ref_val(1, n1+1),
                                                                FEAT->get_ref_val(1, n2+1));
    }
}

double UnivariateNonLinearSRegDiscreteCostFunction::computeUnaryCost(int node, int label) {

    get_target_data(node,estimate_rotation_matrix(_CPgrid.get_coord(node),(*ROTATIONS)[node] * _labels[label]));

    double cost = AbsoluteWeights(node + 1) *
            sim.get_sim_for_min(_sourcedata[node], _targetdata[node],_weights[node]);

    if(_simmeasure == 1 || _simmeasure == 2)
        return cost;
    else
        return -cost;
}

void MultivariateNonLinearSRegDiscreteCostFunction::initialize(int numNodes, int numLabels, int numPairs, int numTriplets) {
    NonLinearSRegDiscreteCostFunction::initialize(numNodes, numLabels, numPairs, numTriplets);
    _sourcedata.clear(); _sourcedata.resize(_SOURCE.nvertices());
    _sourceinrange.clear(); _sourceinrange.resize(_CPgrid.nvertices());
    _targetdata.clear(); _targetdata.resize(_SOURCE.nvertices());
    _weights.clear(); _weights.resize(_SOURCE.nvertices());
}

void MultivariateNonLinearSRegDiscreteCostFunction::get_source_data() {

    for (auto& v: _sourcedata) v.clear();

    for (int i = 0; i < _SOURCE.nvertices(); i++)
    {
        _sourcedata[i].resize(FEAT->get_dim());
        _weights[i].resize(FEAT->get_dim());

        for (int k = 0; k < _CPgrid.nvertices(); k++)
            if (within_controlpt_range(k, i))
                _sourceinrange[k].push_back(i);

        for (int d = 1; d <= FEAT->get_dim(); d++)
        {
            _sourcedata[i][d-1] = FEAT->get_input_val(d, i + 1);
            if (_HIGHREScfweight.Nrows() >= d)
                _weights[i][d-1] = _HIGHREScfweight(d, i + 1);
            else
                _weights[i][d-1] = 1.0;
        }
    }
    resample_weights();
}

void MultivariateNonLinearSRegDiscreteCostFunction::get_target_data(int node, const NEWMAT::Matrix& PtROTATOR) {

    for(unsigned int i = 0; i < _sourceinrange[node].size(); i++)
    {
        _targetdata[_sourceinrange[node][i]].clear();
        _targetdata[_sourceinrange[node][i]].resize(FEAT->get_dim());

        newresampler::Point tmp = PtROTATOR * _SOURCE.get_coord(_sourceinrange[node][i]);

        newresampler::Triangle closest_triangle = targettree->get_closest_triangle(tmp);

        newresampler::Point v0 = closest_triangle.get_vertex_coord(0),
                            v1 = closest_triangle.get_vertex_coord(1),
                            v2 = closest_triangle.get_vertex_coord(2);
                        int n0 = closest_triangle.get_vertex_no(0),
                            n1 = closest_triangle.get_vertex_no(1),
                            n2 = closest_triangle.get_vertex_no(2);

        for(int dim = 0; dim < FEAT->get_dim(); ++dim)
            _targetdata[_sourceinrange[node][i]][dim] = newresampler::barycentric_interpolation(v0, v1, v2, tmp,
                                                                                        FEAT->get_ref_val(dim+1, n0+1),
                                                                                        FEAT->get_ref_val(dim+1, n1+1),
                                                                                        FEAT->get_ref_val(dim+1, n2+1));
    }
}

double MultivariateNonLinearSRegDiscreteCostFunction::computeUnaryCost(int node, int label) {

    double cost	= 0.0;

    get_target_data(node,estimate_rotation_matrix(_CPgrid.get_coord(node),(*ROTATIONS)[node] * _labels[label]));

    for(const auto& i : _sourceinrange[node])
        cost += sim.get_sim_for_min(_sourcedata[i],
                                    _targetdata[i],
                                    _weights[i]);

    if(!_sourceinrange[node].empty()) cost /= (_sourceinrange[node].size());

    if(_simmeasure==1 || _simmeasure==2)
        return AbsoluteWeights(node + 1) * cost;
    else
        return -AbsoluteWeights(node + 1) * cost;
}

void HOUnivariateNonLinearSRegDiscreteCostFunction::initialize(int numNodes, int numLabels, int numPairs, int numTriplets) {
        NonLinearSRegDiscreteCostFunction::initialize(numNodes, numLabels, numPairs, numTriplets);
        _sourcedata.clear(); _sourcedata.resize(_CPgrid.ntriangles());
        _sourceinrange.clear(); _sourceinrange.resize(_CPgrid.ntriangles());
        _targetdata.clear(); _targetdata.resize(_CPgrid.ntriangles());
        _weights.clear(); _weights.resize(_CPgrid.ntriangles());
}

void HOUnivariateNonLinearSRegDiscreteCostFunction::get_source_data() {

    for(auto& v : _sourcedata) v.clear();

    newresampler::Octree cp_tree(_CPgrid);

    for(int i = 0; i < _SOURCE.nvertices(); ++i)
    {
        int closest_triangle = cp_tree.get_closest_triangle(_SOURCE.get_coord(i)).get_no();
        _sourceinrange[closest_triangle].push_back(i);
        _sourcedata[closest_triangle].push_back(FEAT->get_input_val(1, i + 1));
        if(_HIGHREScfweight.Nrows()>=1)
            _weights[closest_triangle].push_back(_HIGHREScfweight(1, i + 1));
        else
            _weights[closest_triangle].push_back(1.0);
    }
    resample_weights();
}

void HOUnivariateNonLinearSRegDiscreteCostFunction::get_target_data(int triplet,
                                                                    const newresampler::Point& new_CP0,
                                                                    const newresampler::Point& new_CP1,
                                                                    const newresampler::Point& new_CP2) {

    _targetdata[triplet].clear();
    _targetdata[triplet].resize(_sourceinrange[triplet].size());

    newresampler::Point CP0 = _CPgrid.get_coord(_triplets[3*triplet]);
    newresampler::Point CP1 = _CPgrid.get_coord(_triplets[3*triplet+1]);
    newresampler::Point CP2 = _CPgrid.get_coord(_triplets[3*triplet+2]);

    #pragma omp parallel for num_threads(mcmc_threads)
    for(int i = 0; i < _sourceinrange[triplet].size(); ++i)
    {
        //newresampler::Point SP = _SOURCE.get_coord(_sourceinrange[triplet][i]);
        newresampler::Point SP = newresampler::project_point(CP0, CP1, CP2, _SOURCE.get_coord(_sourceinrange[triplet][i]));
        newresampler::Point tmp = newresampler::barycentric(CP0,CP1,CP2,SP,new_CP0,new_CP1,new_CP2);
        tmp.normalize();
        tmp *= RAD;

        newresampler::Triangle closest_triangle = targettree->get_closest_triangle(tmp);

        newresampler::Point v0 = closest_triangle.get_vertex_coord(0),
                            v1 = closest_triangle.get_vertex_coord(1),
                            v2 = closest_triangle.get_vertex_coord(2);
                        int n0 = closest_triangle.get_vertex_no(0),
                            n1 = closest_triangle.get_vertex_no(1),
                            n2 = closest_triangle.get_vertex_no(2);

        _targetdata[triplet][i] = newresampler::barycentric_interpolation(v0, v1, v2, tmp,
                                                                   FEAT->get_ref_val(1, n0 + 1),
                                                                   FEAT->get_ref_val(1, n1 + 1),
                                                                   FEAT->get_ref_val(1, n2 + 1));
    }
}

double HOUnivariateNonLinearSRegDiscreteCostFunction::triplet_likelihood(int triplet, int CP_id0, int CP_id1, int CP_id2,
                                                                         const newresampler::Point& CP_def0,
                                                                         const newresampler::Point& CP_def1,
                                                                         const newresampler::Point& CP_def2) {

    get_target_data(triplet, CP_def0, CP_def1, CP_def2);

    double cost = (AbsoluteWeights(CP_id0+1) + AbsoluteWeights(CP_id1+1) + AbsoluteWeights(CP_id2+1)) / 3.0 *
            sim.get_sim_for_min(_sourcedata[triplet], _targetdata[triplet],_weights[triplet]);

    if(_simmeasure == 1 || _simmeasure == 2)
        return cost;
    else
        return -cost;
}

void HOMultivariateNonLinearSRegDiscreteCostFunction::initialize(int numNodes, int numLabels, int numPairs, int numTriplets) {
    NonLinearSRegDiscreteCostFunction::initialize(numNodes, numLabels, numPairs, numTriplets);
    _sourcedata.clear(); _sourcedata.resize(_SOURCE.nvertices());
    _sourceinrange.clear(); _sourceinrange.resize(_CPgrid.ntriangles());
    _targetdata.clear(); _targetdata.resize(_SOURCE.nvertices());
    _weights.clear(); _weights.resize(_SOURCE.nvertices());
}

void HOMultivariateNonLinearSRegDiscreteCostFunction::get_source_data() {

    for(auto& v : _sourcedata) v.clear();

    newresampler::Octree cp_tree(_CPgrid);

    for (int i = 0; i < _SOURCE.nvertices(); ++i)
    {
        _sourceinrange[cp_tree.get_closest_triangle(_SOURCE.get_coord(i)).get_no()].push_back(i);
        _sourcedata[i].resize(FEAT->get_dim());
        _weights[i].resize(FEAT->get_dim());

        for(int dim = 1; dim <= FEAT->get_dim(); ++dim)
        {
            _sourcedata[i][dim-1] = FEAT->get_input_val(dim, i+1);
            if(_HIGHREScfweight.Nrows()>=dim)
                _weights[i][dim-1] = _HIGHREScfweight(dim, i+1);
            else
                _weights[i][dim-1] = 1.0;
        }
    }
    resample_weights();
}

void HOMultivariateNonLinearSRegDiscreteCostFunction::get_target_data(int triplet,
                                                                      const newresampler::Point& new_CP0,
                                                                      const newresampler::Point& new_CP1,
                                                                      const newresampler::Point& new_CP2) {

    newresampler::Point CP0 = _CPgrid.get_coord(_triplets[3*triplet]);
    newresampler::Point CP1 = _CPgrid.get_coord(_triplets[3*triplet+1]);
    newresampler::Point CP2 = _CPgrid.get_coord(_triplets[3*triplet+2]);

    #pragma omp parallel for num_threads(mcmc_threads)
    for(int i = 0; i < _sourceinrange[triplet].size(); ++i)
    {
        _targetdata[_sourceinrange[triplet][i]].clear();
        _targetdata[_sourceinrange[triplet][i]].resize(FEAT->get_dim());

        //newresampler::Point SP = _SOURCE.get_coord(_sourceinrange[triplet][i]);
        newresampler::Point SP = newresampler::project_point(CP0, CP1, CP2, SP);
        newresampler::Point tmp = newresampler::barycentric(CP0,CP1,CP2,SP,new_CP0,new_CP1,new_CP2);
        tmp.normalize();
        tmp *= RAD;

        newresampler::Triangle closest_triangle = targettree->get_closest_triangle(tmp);

        newresampler::Point v0 = closest_triangle.get_vertex_coord(0),
                            v1 = closest_triangle.get_vertex_coord(1),
                            v2 = closest_triangle.get_vertex_coord(2);
                        int n0 = closest_triangle.get_vertex_no(0),
                            n1 = closest_triangle.get_vertex_no(1),
                            n2 = closest_triangle.get_vertex_no(2);

        for(int dim = 0; dim < FEAT->get_dim(); ++dim)
            _targetdata[_sourceinrange[triplet][i]][dim] = newresampler::barycentric_interpolation(v0, v1, v2, tmp,
                                                                         FEAT->get_ref_val(dim + 1, n0 + 1),
                                                                         FEAT->get_ref_val(dim + 1, n1 + 1),
                                                                         FEAT->get_ref_val(dim + 1, n2 + 1));
    }
}

double HOMultivariateNonLinearSRegDiscreteCostFunction::triplet_likelihood(int triplet, int CP_id0, int CP_id1, int CP_id2,
                                                                           const newresampler::Point& CP_def0,
                                                                           const newresampler::Point& CP_def1,
                                                                           const newresampler::Point& CP_def2) {

    get_target_data(triplet, CP_def0, CP_def1, CP_def2);

    double cost = 0.0;

    for(const auto& i : _sourceinrange[triplet])
        cost += sim.get_sim_for_min(_sourcedata[i],
                                    _targetdata[i],
                                    _weights[i]);

    if(_sourceinrange[triplet].size() > 0) cost /= _sourceinrange[triplet].size();

    if(_simmeasure==1 || _simmeasure==2)
        return (AbsoluteWeights(CP_id0+1)+AbsoluteWeights(CP_id1+1)+AbsoluteWeights(CP_id2+1)) / 3.0 * cost;
    else
        return -(AbsoluteWeights(CP_id0+1)+AbsoluteWeights(CP_id1+1)+AbsoluteWeights(CP_id2+1)) / 3.0 * cost;
}

} //namespace newmeshreg
