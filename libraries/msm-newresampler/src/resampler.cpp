/*
Copyright (c) 2022 King's College London, MeTrICS Lab, Renato Besenczi

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

Adaptive Barycentric code used with the permission of Tim Coalson under the same licence as above.
Copyright (C) 2014  Washington University School of Medicine
Original: https://github.com/Washington-University/workbench/blob/master/src/Files/SurfaceResamplingHelper.cxx
*/
#include "resampler.h"

namespace newresampler {

Mesh Resampler::barycentric_data_interpolation(const Mesh& metric_in, const Mesh& sphLow, int nthreads,
                                                                  std::shared_ptr<Mesh> EXCL){

    if (EXCL && EXCL->nvertices() != metric_in.nvertices())
        throw MeshException("Exclusion mask differs in nvertices from data");

    Mesh exclusion = sphLow, interpolated_mesh = sphLow;
    interpolated_mesh.initialize_pvalues(metric_in.get_dimension());
    std::vector<std::map<int,double>> weights = get_adaptive_barycentric_weights(metric_in, sphLow, nthreads, EXCL);

    for(int feat_dim = 0; feat_dim < metric_in.get_dimension(); ++feat_dim) {
#pragma omp parallel for num_threads(nthreads)
        for (int k = 0; k < interpolated_mesh.nvertices(); k++)
        {
            double val = 0.0;

            for (const auto &it: weights[k])
                if (!EXCL || EXCL->get_pvalue(it.first) != 0)
                    val += metric_in.get_pvalue(it.first, feat_dim) * it.second;

            interpolated_mesh.set_pvalue(k, val, feat_dim);
        }
    }

    if (EXCL) {
#pragma omp parallel for num_threads(nthreads)
        for (int k = 0; k < exclusion.nvertices(); k++)
        {
            double excl_val = 0.0;

            for (const auto &it: weights[k])
                if (EXCL->get_pvalue(it.first) != 0)
                    excl_val += EXCL->get_pvalue(it.first) * it.second;

            exclusion.set_pvalue(k, excl_val);
        }
        *EXCL = exclusion;
    }

    return interpolated_mesh;
}

std::vector<std::map<int,double>> Resampler::get_adaptive_barycentric_weights(const Mesh& in_mesh, const Mesh& sphLow, int nthreads, std::shared_ptr<Mesh> EXCL){

    Octree octreeSearch_in(in_mesh);
    std::vector<std::map<int,double>> forward = get_barycentric_weights(sphLow, in_mesh, octreeSearch_in);

    Octree octreeSearch_ref(sphLow);
    std::vector<std::map<int,double>> reverse = get_barycentric_weights(in_mesh, sphLow, octreeSearch_ref);

    std::vector<std::map<int,double>> reverse_reorder, adapt_weights;

    int numOldNodes = in_mesh.nvertices();
    int numNewNodes = sphLow.nvertices();
    std::vector<double> newAreas(numNewNodes, 0.0);
    std::vector<double> oldAreas(numOldNodes, 0.0);
    std::vector<double> correction(numOldNodes, 0.0);

    reverse_reorder.resize(forward.size());
    adapt_weights.resize(forward.size());

    for (int oldNode = 0; oldNode < numOldNodes; ++oldNode)
    {
        oldAreas[oldNode] = compute_vertex_area(oldNode, in_mesh);
        for(auto iter = reverse[oldNode].begin();
            iter != reverse[oldNode].end(); ++iter)
            reverse_reorder[iter->first][oldNode] = iter->second; //this loop can't be parallelized
    }

#pragma omp parallel for num_threads(nthreads)
    for (int newNode = 0; newNode < numNewNodes; ++newNode)
    {
        if(!EXCL || EXCL->get_pvalue(octreeSearch_in.get_closest_vertex_ID(sphLow.get_coord(newNode))) != 0)
        {
            newAreas[newNode] = compute_vertex_area(newNode, sphLow);
            if (reverse_reorder[newNode].size() <= forward[newNode].size())
                // choose whichever has the most samples as this will be the higher res mesh (this differs from source code)
                adapt_weights[newNode] = forward[newNode];
            else
                adapt_weights[newNode] = reverse_reorder[newNode];

            for (auto iter = adapt_weights[newNode].begin();
                 iter != adapt_weights[newNode].end(); ++iter) {
                //begin area correction by multiplying by target node calc_area
                iter->second *= newAreas[newNode];//begin the process of calc_area correction by multiplying by gathering node areas
                correction[iter->first] += iter->second;//now, sum the scattering weights to prepare for first normalization
            }
        }
    }

#pragma omp parallel for num_threads(nthreads)
    for (int newNode = 0; newNode < numNewNodes; ++newNode)
    {
        if(!EXCL || EXCL->get_pvalue(octreeSearch_in.get_closest_vertex_ID(sphLow.get_coord(newNode))) != 0)
        {
            double weight_sum = 0.0;
            for (auto &iter: adapt_weights[newNode])//begin area correction by multiplying by target node calc_area
            {
                iter.second *= oldAreas[iter.first] / correction[iter.first];
                //divide the weights by their scatter sum, then multiply by current areas
                weight_sum += iter.second;//and compute the sum
            }
            if (weight_sum != 0.0)
                //this shouldn't happen unless no nodes remain due to roi, or node areas can be zero
                for (auto &iter: adapt_weights[newNode])
                    iter.second /= weight_sum;//and normalize to a sum of 1
        }
    }

    return adapt_weights;
}

std::vector<std::map<int,double>> Resampler::get_barycentric_weights(const Mesh& low, const Mesh& orig, const Octree& oct, int nthreads){

    std::vector<std::map<int,double>> weights;
    weights.resize(low.nvertices());

#pragma omp parallel for num_threads(nthreads)
    for (int k = 0; k < low.nvertices(); k++)
    {
        Point v0, v1, v2;
        int n0, n1, n2;

        Point ci = low.get_coord(k);
        Triangle closest = oct.get_closest_triangle(ci);

        n0 = closest.get_vertex_no(0);
        n1 = closest.get_vertex_no(1);
        n2 = closest.get_vertex_no(2);
        v0 = closest.get_vertex_coord(0);
        v1 = closest.get_vertex_coord(1);
        v2 = closest.get_vertex_coord(2);

        weights[k] = calc_barycentric_weights(v0, v1, v2, ci, n0, n1, n2);
    }

    return weights;
}

Mesh smooth_data(Mesh& orig, const Mesh& sphLow, double sigma, int nthreads, std::shared_ptr<Mesh> EXCL) {
//---SMOOTHING---//
    check_scale(orig, sphLow);

    NEWMAT::Matrix newdata(orig.get_dimension(), sphLow.nvertices()); newdata = 0;
    Mesh exclusion = sphLow, smoothed = sphLow;
    const double RAD = 100.0;
    const double ang = 4 * asin(sigma / (2 * RAD));

    Octree oct_search(orig);

#pragma omp parallel for num_threads(nthreads)
    for (int i = 0; i < sphLow.nvertices(); i++)
    {
        exclusion.set_pvalue(i,0);
        Point ci = sphLow.get_coord(i);

        int closest_vertex = oct_search.get_closest_vertex_ID(ci);
        Point ref = sphLow.get_coord(closest_vertex);
        ref.normalize();

        double SUM = 0.0, excl_sum = 0.0;

        std::vector<std::pair<int, double>> neighbourhood;

        for(int n = 0; n < sphLow.nvertices(); n++)
        {
            Point actual = sphLow.get_coord(n);
            actual.normalize();
            if ((actual | ref) >= cos(ang))
                neighbourhood.emplace_back(std::make_pair(n, (ref - actual).norm()));
        }

        if(!EXCL || (EXCL->get_pvalue(closest_vertex) > 0))
        {
            for (const auto& neighbour : neighbourhood)
            {
                double geodesic_dist = 2 * RAD * asin(neighbour.second / (2 * RAD));
                double weight = (1 / sqrt(2 * M_PI * sigma * sigma)) * exp(-(geodesic_dist * geodesic_dist) / (2 * sigma * sigma));
                excl_sum += weight;

                if(EXCL) weight = EXCL->get_pvalue(neighbour.first) * weight;

                SUM += weight;
                for(int feat_dim = 0; feat_dim < orig.get_dimension(); ++feat_dim)
                    newdata(feat_dim + 1, i + 1) += orig.get_pvalue(neighbour.first, feat_dim) * weight;
            }

            if(excl_sum != 0.0 && EXCL) exclusion.set_pvalue(i, SUM / excl_sum);

            for(int feat_dim = 1; feat_dim <= orig.get_dimension(); ++feat_dim)
                if(SUM != 0.0) newdata(feat_dim, i + 1) /= SUM;
        }
        else
            exclusion.set_pvalue(i, 0);
    }

    if(EXCL) *EXCL = exclusion;

    smoothed.set_pvalues(newdata);

    return smoothed;
}

Mesh nearest_neighbour_interpolation(Mesh& orig, const Mesh& sphLow, int nthreads, std::shared_ptr<Mesh> EXCL) {
//---NN INTERPOLATION---//
    check_scale(orig, sphLow);

    Mesh exclusion = sphLow, interpolated = sphLow;
    interpolated.initialize_pvalues(orig.get_dimension());

    Octree oct_search(orig);

#pragma omp parallel for num_threads(nthreads)
    for (int i = 0; i < sphLow.nvertices(); i++)
    {
        exclusion.set_pvalue(i, 0);
        int closest_vertex = oct_search.get_closest_vertex_ID(sphLow.get_coord(i));

        if(!EXCL || EXCL->get_pvalue(closest_vertex) != 0)
        {
            if(EXCL) exclusion.set_pvalue(i, EXCL->get_pvalue(closest_vertex));
            for(int feat_dim = 0; feat_dim < orig.get_dimension(); ++feat_dim)
                interpolated.set_pvalue(i, orig.get_pvalue(closest_vertex, feat_dim), feat_dim);
        }
    }

    if(EXCL) *EXCL = exclusion;

    return interpolated;
}

Mesh project_mesh(const Mesh& orig, const Mesh& target, const Mesh& anat, int nthreads) {
//---ANATOMICAL MESH PROJECTION---//
    Resampler resampler(Method::BARY);
    Mesh TRANS = orig;
    Octree octreeSearch(target);

    std::vector<std::map<int,double>> weights = resampler.get_barycentric_weights(orig, target, octreeSearch);

#pragma omp parallel for num_threads(nthreads)
    for (int i = 0; i < orig.nvertices(); i++)
    {
        Point new_coord;
        for (const auto& iter : weights[i])
        {
            if(anat.nvertices() == target.nvertices())
                new_coord += anat.get_coord(iter.first) * iter.second;
            else
                new_coord += target.get_coord(iter.first) * iter.second;
        }
        TRANS.set_coord(i, new_coord);
    }
    return TRANS;
}

Mesh surface_resample(const Mesh& anatOrig, const Mesh& sphOrig, const Mesh& sphLow, int nthreads) {
//---ANATOMICAL MESH RESAMPLING---//
    Resampler resampler(Method::BARY);
    Mesh anatLow = sphLow;
    Octree octreeSearch(sphOrig);

    std::vector<std::map<int,double>> weights = resampler.get_barycentric_weights(sphLow, sphOrig, octreeSearch);

#pragma omp parallel for num_threads(nthreads)
    for(int i = 0; i < sphLow.nvertices(); i++)
    {
        Point newPt;
        for (const auto& it: weights[i])
            newPt += anatOrig.get_coord(it.first) * it.second;
        anatLow.set_coord(i, newPt);
    }

    return anatLow;
}

Mesh metric_resample(const Mesh& metric_in, const Mesh& sphLow, int nthreads, std::shared_ptr<Mesh> EXCL) {
//---METRIC FILE RESAMPLING---//
    Resampler resampler(Method::ADAP_BARY);

    return resampler.barycentric_data_interpolation(metric_in, sphLow, nthreads, EXCL);
}

void barycentric_mesh_interpolation(Mesh& SPH_up, const Mesh& SPH_low_init, const Mesh& SPH_low_final, int nthreads) {

    Resampler R;
    Octree octree_search(SPH_low_init);
    std::vector<std::map<int,double>> weights = R.get_barycentric_weights(SPH_up, SPH_low_init, octree_search);

#pragma omp parallel for num_threads(nthreads)
    for(int i = 0; i < SPH_up.nvertices(); i++)
    {
        Point newPt;
        for(const auto& it : weights[i])
            newPt += SPH_low_final.get_coord(it.first) * it.second;

        newPt.normalize();
        newPt *= 100;
        SPH_up.set_coord(i, newPt);
    }
}

} //namespace newresampler
