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

Adaptive Barycentric code used with the permission of
Tim Coulson under the same licence as above.
Copyright (C) 2014  Washington University School of Medicine
*/
#ifndef NEWRESAMPLER_RESAMPLER_H
#define NEWRESAMPLER_RESAMPLER_H

#include <omp.h>

#include "miscmaths/bfmatrix.h"

#include "mesh.h"
#include "octree.h"

namespace newresampler {

enum class Method { BARY, ADAP_BARY };

class Resampler {

    Method method;

public:
    Resampler(): method(Method::ADAP_BARY) {}
    explicit Resampler(Method m): method(m) {}

    //---ACCESS---//
    void set_method(const Method& m){ method = m; }
    Method get_method() const { return method; };

    //---RESAMPLING AND CALC WEIGHTS, OCTREE-BASED METHODS---//
    Mesh barycentric_data_interpolation(const Mesh& metric_in, const Mesh& sphLow, int nthreads = 1, std::shared_ptr<Mesh> = std::shared_ptr<Mesh>());
    std::vector<std::map<int,double>> get_adaptive_barycentric_weights(const Mesh& in_mesh, const Mesh& sphLow, int nthreads = 1, std::shared_ptr<Mesh> = std::shared_ptr<Mesh>());
    std::vector<std::map<int,double>> get_barycentric_weights(const Mesh& low, const Mesh& orig, const Octree& oct, int nthreads = 1);
};
//---UTILITY---//
Mesh project_mesh(const Mesh& orig, const Mesh& target, const Mesh& anat, int nthreads = 1);

//---ENTRY POINTS FOR RESAMPLER---//
Mesh surface_resample(const Mesh&, const Mesh&, const Mesh&, int nthreads = 1);
Mesh metric_resample(const Mesh&, const Mesh&, int nthreads = 1, std::shared_ptr<Mesh> = std::shared_ptr<Mesh>());
Mesh smooth_data(Mesh& orig, const Mesh& sphLow, double sigma, int nthreads = 1, std::shared_ptr<Mesh> EXCL = std::shared_ptr<Mesh>());
Mesh nearest_neighbour_interpolation(Mesh &orig, const Mesh &sphLow, int nthreads = 1, std::shared_ptr<Mesh> EXCL = std::shared_ptr<Mesh>());
void barycentric_mesh_interpolation(Mesh& SPH_up, const Mesh& SPH_low_init, const Mesh& SPH_low_final, int nthreads = 1);

} //namespace newresampler

#endif //NEWRESAMPLER_RESAMPLER_H
