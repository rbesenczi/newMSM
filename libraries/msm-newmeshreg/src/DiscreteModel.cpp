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
#include "DiscreteModel.h"

namespace newmeshreg {

void NonLinearSRegDiscreteModel::set_parameters(myparam& PAR) {
    myparam::iterator it;
    it=PAR.find("dOPT"); optimiser = std::get<std::string>(it->second);
    it=PAR.find("SGres"); m_SGres = std::get<int>(it->second);
    it=PAR.find("regularisermode"); m_regoption = std::get<int>(it->second);
    it=PAR.find("multivariate"); m_multivariate = std::get<bool>(it->second);
    it=PAR.find("verbosity"); m_verbosity = std::get<bool>(it->second);
    it=PAR.find("outdir"); m_outdir = std::get<std::string>(it->second);
    it=PAR.find("TriLikelihood"); m_triclique = std::get<bool>(it->second);
    it=PAR.find("rescalelabels"); m_rescalelabels = std::get<bool>(it->second);
    it=PAR.find("numthreads"); _nthreads = std::get<int>(it->second);
    it=PAR.find("range"); range=std::get<double>(it->second);
    if(m_regoption == 1) _pairwise = true;
}

void NonLinearSRegDiscreteModel::initialize_cost_function(bool MV, myparam& P) {

    if (MV)
        if(m_triclique)
            costfct = std::shared_ptr<NonLinearSRegDiscreteCostFunction>(new HOMultivariateNonLinearSRegDiscreteCostFunction());
        else
            costfct = std::shared_ptr<NonLinearSRegDiscreteCostFunction>(new MultivariateNonLinearSRegDiscreteCostFunction());
    else
        if (m_triclique)
            costfct = std::shared_ptr<NonLinearSRegDiscreteCostFunction>(new HOUnivariateNonLinearSRegDiscreteCostFunction());
        else
            costfct = std::shared_ptr<NonLinearSRegDiscreteCostFunction>(new UnivariateNonLinearSRegDiscreteCostFunction());

    costfct->set_parameters(P);
}

void NonLinearSRegDiscreteModel::Initialize(const newresampler::Mesh& CONTROLGRID) {

    m_CPgrid = CONTROLGRID;

    //---SET LOW RES DEFORMATION GRID & INITIALISE ASSOCIATED MRF PARAMS---//
    m_num_nodes = m_CPgrid.nvertices();
    m_num_pairs = 0;
    initLabeling();

    //---CALCULATE (MEAN) VERTEX SEPARATIONS FOR EACH VERTEX
    NEWMAT::ColumnVector vMAXmvd;
    vMAXmvd.ReSize(m_CPgrid.nvertices());
    vMAXmvd = 0;

    #pragma omp parallel for num_threads(_nthreads)
    for (int k = 0; k < m_CPgrid.nvertices(); k++)
    {
        newresampler::Point CP = m_CPgrid.get_coord(k);
        for (auto it = m_CPgrid.nbegin(k); it != m_CPgrid.nend(k); it++) {
            double dist = 2 * RAD * asin((CP-m_CPgrid.get_coord(*it)).norm() / (2 * RAD));
            if(dist > vMAXmvd(k+1)) vMAXmvd(k+1) = dist;
        }
    }

    double MVDmax = m_CPgrid.calculate_MaxVD();

    m_maxs_dist = _labeldist * MVDmax;

    //---INITIALIZE COSTFCT---//
    costfct->set_meshes(m_TARGET, m_SOURCE, m_CPgrid, 0);
    costfct->set_spacings(vMAXmvd, MVDmax);

    m_iter = 1;
    m_scale = 1.0;

    if (_pairwise)
        estimate_pairs();
    else
        estimate_triplets();

    //---INITIALIAZE LABEL GRID---//
    Initialize_sampling_grid();

    //---INITIALIZE NEIGHBOURHOODS---//
    m_inputtree = std::make_shared<newresampler::Octree>(m_TARGET);
}

void NonLinearSRegDiscreteModel::Initialize_sampling_grid() {
    //---LABELS USING HIGHER RES GRID---//
    m_samplinggrid = newresampler::make_mesh_from_icosa(m_SGres);
    true_rescale(m_samplinggrid,RAD);
    // find the first centroid with 6 neighbours
    for (int i = 0; i < m_samplinggrid.nvertices(); i++)
        if (m_samplinggrid.get_total_neighbours(i) == 6)
        {
            m_centroid = i;
            break;
        }
    label_sampling_grid(m_centroid,m_maxs_dist,m_samplinggrid);
}

void NonLinearSRegDiscreteModel::label_sampling_grid(int centroid, double dist, newresampler::Mesh& Grid) {

    m_samples.clear();
    m_barycentres.clear();
    std::vector<int> getneighbours, newneighbours;
    NEWMAT::ColumnVector found(Grid.nvertices()), found_tr(Grid.ntriangles());
    found = 0; found_tr = 0;

    centre = Grid.get_coord(centroid);
    getneighbours.push_back(centroid);
    m_samples.push_back(centre);
    m_barycentres.push_back(centre);

    // searches for neighbours of the centroid that are within the max sampling distance
    while(!getneighbours.empty())
    {
        for(int& getneighbour : getneighbours)
        {
            for(auto j = Grid.nbegin(getneighbour); j != Grid.nend(getneighbour); j++)
            {
                newresampler::Point sample = Grid.get_coord(*j);
                if((sample-centre).norm() <= dist && (found(*j+1)==0) && *j != centroid)
                {
                    m_samples.push_back(Grid.get_coord(*j));  // pt-centroid equals deformation vector
                    newneighbours.push_back(*j);
                    found(*j+1) = 1;
                }
            }

            for(auto j = Grid.tIDbegin(getneighbour); j != Grid.tIDend(getneighbour); j++)
            {
                newresampler::Point v1 = Grid.get_triangle_vertex(*j,0),
                                    v2 = Grid.get_triangle_vertex(*j,1),
                                    v3 = Grid.get_triangle_vertex(*j,2);
                newresampler::Point bary((v1.X + v2.X + v3.X) / 3,
                                         (v1.Y + v2.Y + v3.Y) / 3,
                                         (v1.Z + v2.Z + v3.Z) / 3);
                bary.normalize();
                bary = bary * RAD;

                if((bary - centre).norm() <= dist && (bary - centre).norm() > 0 && found_tr(*j+1) == 0)
                {
                    for (auto& m_barycentre: m_barycentres)
                        if (abs(1 - (((bary - centre) | (m_barycentre - centre)) /
                                     ((bary - centre).norm() * (m_barycentre - centre).norm()))) < 1e-2)
                            found_tr(*j + 1) = 1;

                    if (found_tr(*j + 1) == 0)
                        m_barycentres.push_back(bary);

                    found_tr(*j+1) = 1;
                }
            }
        }
        getneighbours = newneighbours;
        newneighbours.clear();
    }
}

std::vector<newresampler::Point> NonLinearSRegDiscreteModel::rescale_sampling_grid() {

    std::vector<newresampler::Point> newlabels(m_samples.size());

    if(m_verbosity) std::cout << " resample labels " << m_scale << " length scale " << (centre-m_samples[1]).norm() << std::endl;

    if(m_scale >= 0.25)
    {
        for (int i = 0; i < m_samples.size(); i++)
        {
            newresampler::Point newsample = centre + (centre - m_samples[i]) * m_scale;
            newsample.normalize();
            newlabels[i] = newsample * 100;
        }
    }
    else
    {
        m_scale = 1;
        newlabels = m_samples;
    }
    m_scale *= 0.8;
    return newlabels;
}

void NonLinearSRegDiscreteModel::setupCostFunction() {

    resetLabeling(); // initialise label array to zero
    //---use geodesic distances---//
    costfct->reset_CPgrid(m_CPgrid, 0);

    if(m_iter == 1)
    {
        costfct->initialize_regulariser();
        costfct->set_octrees(m_inputtree);
    }

    costfct->reset_anatomical();
    get_rotations();
    // instead of recalulating the source->CP neighbours, these are now constant
    // (as source is moving with CPgrid) we we just need to recalculate
    // the rotations of the label grid to the cp vertices

    if(m_debug)
    {
        m_SOURCE.save(m_outdir + "SOURCE-" + std::to_string(m_iter) + ".surf");
        m_SOURCE.save(m_outdir + "SOURCE-" + std::to_string(m_iter) + ".func");
        if (m_iter == 1) m_TARGET.save(m_outdir + "TARGET.surf");
        m_CPgrid.save(m_outdir + "CPgrid-" + std::to_string(m_iter) + ".surf");
    }

    //---set samples (labels vary from iter to iter)---//
    if (m_rescalelabels)
        m_labels = rescale_sampling_grid();
    else if (m_iter % 2 == 0)
        m_labels = m_samples;
    else
        m_labels = m_barycentres;

    m_num_labels = m_labels.size();

    costfct->set_labels(m_labels,m_ROT);
    if(m_verbosity) std::cout << " initialise cost function " << m_iter <<  std::endl;

    costfct->initialize(m_num_nodes,m_num_labels,m_num_pairs,m_num_triplets);
    costfct->get_source_data();

    if (_pairwise) costfct->setPairs(pairs);
    else costfct->setTriplets(triplets);

    if(optimiser == "MCMC") costfct->set_mcmc_threads(_nthreads);

    m_iter++;
}

void NonLinearSRegDiscreteModel::applyLabeling() {

    // rotate sampling points to overlap with control point transform grid point to new position given by label
    for (int i = 0; i < m_CPgrid.nvertices(); i++)
        m_CPgrid.set_coord(i, m_ROT[i] * m_labels[labeling[i]]);
}

void NonLinearSRegDiscreteModel::estimate_pairs() {

    int pair = 0;

    for (int i = 0; i < m_CPgrid.nvertices(); i++) // estimate the total number of edge pairs
        m_num_pairs += m_CPgrid.get_total_neighbours(i);

    m_num_pairs /= 2;
    pairs = new int[m_num_pairs*2];

    for (int i = 0; i < m_CPgrid.nvertices(); i++)
        for (auto j = m_CPgrid.nbegin(i); j != m_CPgrid.nend(i); j++)
            if (*j > i)
            {
                int node_ids[2] = {i, *j};
                std::sort(std::begin(node_ids), std::end(node_ids));
                pairs[2*pair  ] = node_ids[0];
                pairs[2*pair+1] = node_ids[1];
                pair++;
            }
}

void NonLinearSRegDiscreteModel::estimate_triplets() {

    m_num_triplets = m_CPgrid.ntriangles();
    triplets = new int[m_num_triplets*3];

    for(int i = 0; i < m_CPgrid.ntriangles(); i++)
    {
        int node_ids[3] = {m_CPgrid.get_triangle(i).get_vertex_no(0),
                           m_CPgrid.get_triangle(i).get_vertex_no(1),
                           m_CPgrid.get_triangle(i).get_vertex_no(2) };
        std::sort(std::begin(node_ids), std::end(node_ids));
        triplets[3*i  ] = node_ids[0];
        triplets[3*i+1] = node_ids[1];
        triplets[3*i+2] = node_ids[2];
    }
}

void NonLinearSRegDiscreteModel::get_rotations() {

    // rotates sampling grid to each control point
    m_ROT.clear();
    m_ROT.resize(m_CPgrid.nvertices());

    #pragma omp parallel for num_threads(_nthreads)
    for (int k = 0; k < m_CPgrid.nvertices(); k++)
        m_ROT[k] = estimate_rotation_matrix(centre, m_CPgrid.get_coord(k));
}

} //namespace newmeshreg
