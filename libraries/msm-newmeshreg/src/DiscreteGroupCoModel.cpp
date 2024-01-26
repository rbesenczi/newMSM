#include "DiscreteGroupCoModel.h"

namespace newmeshreg {

void DiscreteGroupCoModel::initialize_pairs() {

    delete[] pairs;
    m_num_pairs = 0;
    for (int subject_A = 0; subject_A < m_num_subjects; subject_A++)
        for (int vertex = 0; vertex < control_grid_size; vertex++)
            for (int subject_B = subject_A + 1; subject_B < m_num_subjects; subject_B++)
                m_num_pairs++;
    pairs = new int[2 * m_num_pairs];
}

void DiscreteGroupCoModel::estimate_pairs() {

    int pair = 0;
    std::vector<std::shared_ptr<newresampler::Octree>> cp_grid_trees(m_num_subjects);

    #pragma omp parallel for num_threads(_nthreads)
    for (int subject = 0; subject < m_num_subjects; ++subject)
        cp_grid_trees[subject] = std::make_shared<newresampler::Octree>(m_controlmeshes[subject]);

    for (int subject_A = 0; subject_A < m_num_subjects; subject_A++)
        for (int vertex = 0; vertex < control_grid_size; vertex++) {
            newresampler::Point CP = m_controlmeshes[subject_A].get_coord(vertex);
            for (int subject_B = subject_A + 1; subject_B < m_num_subjects; subject_B++) {
                pairs[2 * pair] = subject_A * control_grid_size + vertex;
                pairs[2 * pair + 1] =
                        subject_B * control_grid_size + cp_grid_trees[subject_B]->get_closest_vertex_ID(CP);
                pair++;
            }
        }
}

void DiscreteGroupCoModel::estimate_triplets() {

    delete[] triplets;
    triplets = new int[3 * m_num_triplets];

    #pragma omp parallel for num_threads(_nthreads)
    for (int subject = 0; subject < m_num_subjects; subject++)
        for (int triangle = 0; triangle < cp_triangles; triangle++) {
            int node_ids[3] =
                    {m_controlmeshes[subject].get_triangle_vertexID(triangle, 0) + subject * control_grid_size,
                     m_controlmeshes[subject].get_triangle_vertexID(triangle, 1) + subject * control_grid_size,
                     m_controlmeshes[subject].get_triangle_vertexID(triangle, 2) + subject * control_grid_size};
            std::sort(std::begin(node_ids), std::end(node_ids));
            triplets[subject * cp_triangles * 3 + triangle * 3] = node_ids[0];
            triplets[subject * cp_triangles * 3 + triangle * 3 + 1] = node_ids[1];
            triplets[subject * cp_triangles * 3 + triangle * 3 + 2] = node_ids[2];
        }
}

void DiscreteGroupCoModel::get_rotations() {

    m_ROT.clear();
    m_ROT.resize(m_num_subjects * control_grid_size);

    #pragma omp parallel for num_threads(_nthreads)
    for (int subject = 0; subject < m_num_subjects; subject++)
        for (int vertex = 0; vertex < control_grid_size; vertex++)
            m_ROT[subject * control_grid_size + vertex] =
                    newresampler::estimate_rotation_matrix(centre,m_controlmeshes[subject].get_coord(vertex));

    warp_rotations.clear();
    warp_rotations.resize(2);
    warp_rotations[0].resize(warps[0].size());
    warp_rotations[1].resize(warps[1].size());

    for (int group = 0; group < 2; group++)
        for (int subject = 0; subject < warps[group].size(); subject++) {
            warp_rotations[group][subject].resize(control_grid_size);
            for (int vertex = 0; vertex < control_grid_size; vertex++)
                warp_rotations[group][subject][vertex] =
                        newresampler::estimate_rotation_matrix(centre,warps[group][subject].get_coord(vertex));
        }
}

void DiscreteGroupCoModel::get_patch_data() {

    std::vector<std::map<int, float>> patch_data(control_grid_size * m_num_subjects * m_num_labels);

    #pragma omp parallel for num_threads(_nthreads)
    for (int subject = 0; subject < m_num_subjects; subject++)
        for (int label = 0; label < m_num_labels; label++) {
            newresampler::Mesh rotated_mesh = m_datameshes[subject];
            rotated_mesh.set_pvalues(FEAT->get_data_matrix(subject));

            if (label > 0)
                for (int datapoint = 0; datapoint < rotated_mesh.nvertices(); datapoint++)
                    rotated_mesh.set_coord(datapoint,
                                           estimate_rotation_matrix(centre, rotated_mesh.get_coord(datapoint)) *
                                           m_labels[label]); // rigid rotation, bugs expected

            rotated_mesh = newresampler::metric_resample(rotated_mesh, m_template);

            for (int vertex = 0; vertex < control_grid_size; vertex++) {
                std::map<int, float> patchdata;
                const newresampler::Point rotated_CP =
                        m_ROT[subject * control_grid_size + vertex] * m_labels[label];
                for (int datapoint = 0; datapoint < rotated_mesh.nvertices(); datapoint++)
                    if (((2 * RAD * asin((rotated_CP - rotated_mesh.get_coord(datapoint)).norm() / (2 * RAD)))
                         < range * spacings[subject](vertex + 1)))
                        patchdata[datapoint] = rotated_mesh.get_pvalue(datapoint);

                patch_data[subject * control_grid_size * m_num_labels + vertex * m_num_labels + label] = patchdata;
            }
        }
    costfct->set_patch_data(patch_data);
}

void DiscreteGroupCoModel::get_spacings() {

    spacings.clear();
    spacings.resize(m_num_subjects);

    #pragma omp parallel for num_threads(_nthreads)
    for (int subject = 0; subject < m_num_subjects; subject++) {
        NEWMAT::ColumnVector vMAXmvd(control_grid_size);
        vMAXmvd = 0;
        for (int vertex = 0; vertex < control_grid_size; vertex++) {
            newresampler::Point CP = m_controlmeshes[subject].get_coord(vertex);
            for (auto it = m_controlmeshes[subject].nbegin(vertex);
                 it != m_controlmeshes[subject].nend(vertex); it++) {
                double dist = 2 * RAD * asin((CP - m_controlmeshes[subject].get_coord(*it)).norm() / (2 * RAD));
                if (dist > vMAXmvd(vertex + 1)) vMAXmvd(vertex + 1) = dist;
            }
        }
        spacings[subject] = vMAXmvd;
    }
}

void DiscreteGroupCoModel::Initialize(const newresampler::Mesh &controlgrid) {

    m_controlmeshes.clear();
    m_controlmeshes.resize(m_num_subjects, controlgrid);
    control_grid_size = controlgrid.nvertices();
    cp_triangles = controlgrid.ntriangles();
    m_maxs_dist = _labeldist * controlgrid.calculate_MaxVD();
    m_num_nodes = control_grid_size * m_num_subjects;
    m_num_triplets = m_num_subjects * cp_triangles;

    m_iter = 1;

    initLabeling();
    Initialize_sampling_grid();
    initialize_pairs();
    estimate_triplets();

    costfct->set_meshes(m_datameshes, controlgrid, m_num_subjects);
}

void DiscreteGroupCoModel::setupCostFunction() {

    if (m_verbosity) std::cout << "Initialising cost function " << m_iter << std::endl;

    resetLabeling();

    for (int subject = 0; subject < m_num_subjects; subject++)
        costfct->reset_CPgrid(m_controlmeshes[subject], subject);

    estimate_pairs();

    get_spacings();
    get_rotations();

    costfct->set_iter(m_iter);

    if (m_iter % 2 == 0) m_labels = m_samples;
    else m_labels = m_barycentres;
    m_num_labels = (int) m_labels.size();

    costfct->set_labels(m_labels, m_ROT);
    costfct->set_warp_rotations(warp_rotations);
    costfct->set_warps(warps);

    get_patch_data();

    costfct->setPairs(pairs);
    costfct->setTriplets(triplets);

    costfct->initialize(m_num_nodes, m_num_labels, m_num_pairs, m_num_triplets);

    if (m_verbosity)
        std::cout << " numpoints " << m_num_nodes << " m_num_labels " << m_num_labels << " m_num_pairs "
                  << m_num_pairs << std::endl;

    m_iter++;
}

} //namespace meshreg
