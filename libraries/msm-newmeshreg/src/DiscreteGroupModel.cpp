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
#include "DiscreteGroupModel.h"

namespace newmeshreg {

void DiscreteGroupModel::initialize_pairs() {

    delete[] pairs;
    m_num_pairs = 0;
    for (int n = 0; n < m_num_subjects; n++)
        for (int i = 0; i < control_grid_size; i++)
            for (int n2 = n+1; n2 < m_num_subjects; n2++)
                m_num_pairs++;
    pairs = new int[2 * m_num_pairs];
}

void DiscreteGroupModel::estimate_pairs() {

    int pair = 0;
    std::vector<std::shared_ptr<newresampler::Octree>> cp_grid_trees(m_num_subjects);

    #pragma omp parallel for num_threads(_nthreads)
    for (int subject = 0; subject < m_num_subjects; ++subject)
        cp_grid_trees[subject] = std::make_shared<newresampler::Octree>(m_controlmeshes[subject]);

    for (int subject_A = 0; subject_A < m_num_subjects; subject_A++)
        for (int vertex = 0; vertex < control_grid_size; vertex++) {
            newresampler::Point CP = m_controlmeshes[subject_A].get_coord(vertex);
            for (int subject_B = subject_A + 1; subject_B < m_num_subjects; subject_B++) {
                pairs[2*pair  ] = subject_A * control_grid_size + vertex;
                pairs[2*pair+1] = subject_B * control_grid_size + cp_grid_trees[subject_B]->get_closest_vertex_ID(CP);
                pair++;
            }
        }
}

void DiscreteGroupModel::estimate_triplets() {

    int num_triangles = m_controlmeshes[0].ntriangles();
    m_num_triplets = m_num_subjects * num_triangles;
    constexpr int vertex_per_triangle = 3;

    delete[] triplets;
    triplets = new int[vertex_per_triangle * m_num_triplets];

    #pragma omp parallel for num_threads(_nthreads)
    for (int subject = 0; subject < m_num_subjects; subject++)
        for (int triangle = 0; triangle < num_triangles; triangle++)
        {
            int node_ids[vertex_per_triangle] =
                  {m_controlmeshes[subject].get_triangle_vertexID(triangle, 0) + subject * control_grid_size,
                   m_controlmeshes[subject].get_triangle_vertexID(triangle, 1) + subject * control_grid_size,
                   m_controlmeshes[subject].get_triangle_vertexID(triangle, 2) + subject * control_grid_size };
            std::sort(std::begin(node_ids), std::end(node_ids));
            triplets[subject * num_triangles * vertex_per_triangle + triangle * vertex_per_triangle    ] = node_ids[0];
            triplets[subject * num_triangles * vertex_per_triangle + triangle * vertex_per_triangle + 1] = node_ids[1];
            triplets[subject * num_triangles * vertex_per_triangle + triangle * vertex_per_triangle + 2] = node_ids[2];
        }
}

void DiscreteGroupModel::get_rotations(std::vector<NEWMAT::Matrix>& ROT) {

    ROT.clear();
    ROT.resize(m_num_subjects * control_grid_size);
    const newresampler::Point ci = m_samplinggrid.get_coord(m_centroid);

    #pragma omp parallel for num_threads(_nthreads)
    for (int subject = 0; subject < m_num_subjects; subject++)
        for (int vertex = 0; vertex < control_grid_size; vertex++)
            ROT[subject * control_grid_size + vertex] = estimate_rotation_matrix(ci, m_controlmeshes[subject].get_coord(vertex));
}

void DiscreteGroupModel::get_rotated_meshes() {

    rotated_meshes.clear();
    rotated_meshes.resize(m_num_subjects * m_num_labels);

    #pragma omp parallel for num_threads(_nthreads)
    for(int subject = 0; subject < m_num_subjects; subject++) {
        for (int label = 0; label < m_num_labels; label++) {
            newresampler::Mesh rotated_mesh = m_datameshes[subject];
            rotated_mesh.set_pvalues(FEAT->get_data_matrix(subject));
            for (int vertex = 0; vertex < control_grid_size; vertex++) {
                const newresampler::Point CP = m_controlmeshes[subject].get_coord(vertex);
                const NEWMAT::Matrix CP_label_rotation = estimate_rotation_matrix(CP, m_ROT[subject * control_grid_size + vertex] * m_labels[label]);
                for (int datapoint = 0; datapoint < m_datameshes[subject].nvertices(); datapoint++) {
                    const newresampler::Point SP = m_datameshes[subject].get_coord(datapoint);
                    if (((2 * RAD * asin((CP - SP).norm() / (2 * RAD))) < range * spacings[subject](vertex + 1)))
                        rotated_mesh.set_coord(datapoint, SP * CP_label_rotation);
                }
            }
            rotated_meshes[subject*m_num_labels+label] = newresampler::metric_resample(rotated_mesh, m_template).get_pvalues();
        }
    }
}

void DiscreteGroupModel::get_patch_data() {

    patch_data.clear();
    patch_data.resize(control_grid_size * m_num_subjects * m_num_labels);

    #pragma omp parallel for num_threads(_nthreads)
    for (int subject = 0; subject < m_num_subjects; subject++) {
        for (int label = 0; label < m_num_labels; label++) {
            for (int vertex = 0; vertex < control_grid_size; vertex++) {
                std::map<int, double> patchdata;
                const newresampler::Point CP = m_ROT[subject * control_grid_size + vertex] * m_labels[label];
                for (int datapoint = 0; datapoint < m_template.nvertices(); datapoint++)
                    if (((2 * RAD * asin((CP - m_template.get_coord(datapoint)).norm() / (2 * RAD))) < range * spacings[subject](vertex + 1)))
                        patchdata[datapoint] = rotated_meshes[subject * m_num_labels + label](1, datapoint + 1);

                patch_data[subject * control_grid_size * m_num_labels + vertex * m_num_labels + label] = patchdata;
            }
        }
    }
}

void DiscreteGroupModel::Initialize(const newresampler::Mesh& controlgrid) {

    m_controlmeshes.clear();
    m_controlmeshes.resize(m_num_subjects, controlgrid);
    control_grid_size = controlgrid.nvertices();
    m_num_nodes = control_grid_size * m_num_subjects;

    initLabeling();

    spacings.clear();
    spacings.resize(m_num_subjects);
    #pragma omp parallel for num_threads(_nthreads)
    for(int subject = 0; subject < m_num_subjects; subject++)
    {
        NEWMAT::ColumnVector vMAXmvd(control_grid_size);
        vMAXmvd = 0;
        for (int vertex = 0; vertex < control_grid_size; vertex++)
        {
            newresampler::Point CP = m_controlmeshes[subject].get_coord(vertex);
            for (auto it = m_controlmeshes[subject].nbegin(vertex); it != m_controlmeshes[subject].nend(vertex); it++)
            {
                double dist = 2*RAD * asin((CP - m_controlmeshes[subject].get_coord(*it)).norm() / (2 * RAD));
                if(dist > vMAXmvd(vertex + 1)) vMAXmvd(vertex + 1) = dist;
            }
        }
        spacings[subject] = vMAXmvd;
    }

    MVD = controlgrid.calculate_MeanVD();
    m_maxs_dist = _labeldist * controlgrid.calculate_MaxVD();

    costfct->set_meshes(m_template, m_datameshes[0], controlgrid, m_num_subjects);
    //costfct->set_group_spacings(spacings);

    m_iter = 1;

    //---GET BETWEEN MESH GROUPINGS---//
    initialize_pairs();
    estimate_pairs();
    estimate_triplets();

    Initialize_sampling_grid();
    get_rotations(m_ROT);
}

void DiscreteGroupModel::setupCostFunction() {

    //---INIT---//
    if(m_verbosity)
        std::cout << " initialize cost function " << m_iter << " m_num_triplets " << m_num_triplets << std::endl;

    resetLabeling(); // initialise label array to zero

    for (int subject = 0; subject < m_num_subjects; subject++)
        costfct->reset_CPgrid(m_controlmeshes[subject], subject);

    get_rotations(m_ROT);

    costfct->set_iter(m_iter);

    if (m_iter % 2 == 0)
        m_labels = m_samples;
    else
        m_labels = m_barycentres;

    m_num_labels = m_labels.size();

    costfct->set_labels(m_labels,m_ROT);
    costfct->initialize(m_num_nodes, m_num_labels, m_num_pairs, m_num_triplets);

    estimate_pairs();
    get_rotated_meshes();
    get_patch_data();

    costfct->setPairs(pairs);
    costfct->setTriplets(triplets);
    costfct->set_patch_data(patch_data);

    if(m_verbosity)
        std::cout << " numpoints " << m_num_nodes << " m_num_labels " << m_num_labels << " m_num_pairs " << m_num_pairs << std::endl;

    m_iter++;
}

} //namespace newmeshreg
