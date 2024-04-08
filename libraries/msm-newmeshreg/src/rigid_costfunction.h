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
#ifndef NEWMESHREG_RIGID_COSTFUNCTION_H
#define NEWMESHREG_RIGID_COSTFUNCTION_H

#include "featurespace.h"
#include "similarities.h"
#include "newresampler/octree.h"
#include <omp.h>

namespace newmeshreg {

class Rigid_cost_function {
// basic mesh cost function - uses the gradient of the similarity function, estimated using weighted least squares

    newresampler::Mesh TARGET, SOURCE;

    std::shared_ptr<Neighbourhood> nbh;
    std::shared_ptr<newresampler::Octree> targettree;
    std::shared_ptr<featurespace> FEAT; // holds data
    sparsesimkernel sim; // similarity matrix
    double MVD = 0.0, min_sigma = 0.0;
    NEWMAT::ColumnVector current_sim;

    // user defined parameters
    int simmeasure = 1, iters = 20;
    float stepsize = 0.01, spacing = 0.5;
    bool verbosity = false;
    int numthreads = 1;

    //---SIMILARITY GRADIENT ESTIMATION---//
    void WLS_simgradient(const newresampler::Tangs& tangs, int index, const std::vector<int>& querypoints);
    // prepares data for sim gradient calculation
    void Evaluate_SIMGradient(int i, const newresampler::Tangs& tangs);

    //---TRANSFORM AND EVALUATE---//
    void rotate_in_mesh(double a1, double a2, double a3);
    double rigid_cost_mesh(double dw1, double dw2, double dw3);
    bool get_all_neighbours(int index, std::vector<int>& N, int n, const newresampler::Mesh& REF, std::shared_ptr<Neighbourhood>& _rel, MISCMATHS::SpMat<int>& found);

public:
    Rigid_cost_function(newresampler::Mesh target, newresampler::Mesh source,
                        std::shared_ptr<featurespace>& features);
    void set_parameters(myparam& PAR);

    //---INITIALIZE AND UPDATE---//
    void initialise();
    void set_simmeasure(int simval) { simmeasure = simval; }
    void update_source(const newresampler::Mesh& M) { SOURCE = M; }

    //---EXECUTE---//
    newresampler::Mesh run();
};

} //namespace newmeshreg

#endif //NEWMESHREG_RIGID_COSTFUNCTION_H
