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
#include "rigid_costfunction.h"

namespace newmeshreg {

Rigid_cost_function::Rigid_cost_function(
        newresampler::Mesh target,
        newresampler::Mesh source,
        std::shared_ptr<featurespace>& features)
        :TARGET(std::move(target)), SOURCE(std::move(source)), FEAT(features) {}

void Rigid_cost_function::initialise() {

    current_sim.ReSize(SOURCE.nvertices());
    min_sigma = MVD = SOURCE.calculate_MeanVD();
    sim.resize(TARGET.nvertices(), SOURCE.nvertices());  // creates sparse sim kernel to fill with similarities of nearest neighbours only
    sim.set_reference(FEAT->get_reference_data());
    sim.set_input(FEAT->get_input_data());
    sim.initialise(simmeasure);
    nbh = std::make_shared<Neighbourhood>();
    nbh->update(SOURCE, TARGET, 2 * asin(4 * MVD / (2 * RAD)), numthreads);
    sim.set_neighbourhood(nbh);
    targettree = std::make_shared<newresampler::Octree>(TARGET);

    #pragma omp parallel for num_threads(numthreads)
    for (int i = 0; i < SOURCE.nvertices(); i++)
        sim.calculate_sim_column_nbh(i);
}

void Rigid_cost_function::set_parameters(myparam& PAR){
    myparam::iterator it;
    it = PAR.find("iters"); iters = std::get<int>(it->second);
    it = PAR.find("simmeasure"); simmeasure = std::get<int>(it->second);
    it = PAR.find("verbosity"); verbosity = std::get<bool>(it->second);
    it = PAR.find("stepsize"); stepsize = std::get<float>(it->second);
    it = PAR.find("gradsampling"); spacing = std::get<float>(it->second);
    it = PAR.find("numthreads"); numthreads = std::get<int>(it->second);
}

void Rigid_cost_function::WLS_simgradient(const newresampler::Tangs& tangs, int index, const std::vector<int>& querypoints) {

    double x11, x21, y11, y21, SUM = 0.0, JPsim = 0.0;
    newresampler::Point origin = tangs.e1 * tangs.e2;
    origin.normalize(); origin *= RAD;
    newresampler::project_point(SOURCE.get_coord(index) - origin, tangs, y11, y21);

    for(const int& querypoint : querypoints)
    {
        newresampler::Point cr = TARGET.get_coord(querypoint);
        newresampler::project_point(cr - origin, tangs, x11, x21);
        double dist_1 = x11 - y11;
        double dist_2 = x21 - y21;

        if((dist_1 * dist_1 + dist_2 * dist_2) > 0)
        {
            double weight = exp(-(dist_1 * dist_1 + dist_2 * dist_2) / (2 * min_sigma * min_sigma));
            SUM += weight;
            JPsim += sim.peek(querypoint+1,index+1) * weight;
        }
    }

    if(SUM > 0) JPsim /= SUM;

    current_sim(index+1) = JPsim;
}

void Rigid_cost_function::Evaluate_SIMGradient(int i, const newresampler::Tangs& tangs) {

    MISCMATHS::SpMat<int> found(TARGET.nvertices(),1);
    newresampler::Point point = SOURCE.get_coord(i);
    bool update = false;

    if(nbh->nrows(i) > 0)
    {
        std::vector<int> querypoints;
        newresampler::Triangle closest_triangle = targettree->get_closest_triangle(point);

        if(get_all_neighbours(i, querypoints, closest_triangle.get_vertex_no(0), TARGET, nbh, found)) update = true;
        if(get_all_neighbours(i, querypoints, closest_triangle.get_vertex_no(1), TARGET, nbh, found) || update) update = true;
        if(get_all_neighbours(i, querypoints, closest_triangle.get_vertex_no(2), TARGET, nbh, found) || update) update = true;
        if(update)
        {
            nbh->at(i) = querypoints;
            sim.calculate_sim_column_nbh(i);
        }
        WLS_simgradient(tangs,i, querypoints);
    }
}

void Rigid_cost_function::rotate_in_mesh(double a1, double a2, double a3){

    #pragma omp parallel for num_threads(numthreads)
    for (int index = 0; index < SOURCE.nvertices(); index++)
    {
        newresampler::Point cii = SOURCE.get_coord(index);
        NEWMAT::ColumnVector V(3);
        V(1) = cii.X; V(2) = cii.Y; V(3) = cii.Z;
        NEWMAT::ColumnVector VR = newresampler::euler_rotate(V, a1, a2, a3);
        SOURCE.set_coord(index,newresampler::Point(VR(1),VR(2),VR(3)));
    }
}

double Rigid_cost_function::rigid_cost_mesh(double dw1, double dw2, double dw3){

    double SUM = 0.0;
    newresampler::Mesh tmp = SOURCE;

    rotate_in_mesh(dw1, dw2, dw3);

    #pragma omp parallel for num_threads(numthreads)
    for (int index = 0; index < SOURCE.nvertices(); index++)
        Evaluate_SIMGradient(index, calculate_tangs(index, SOURCE));

    for(int i = 1; i <= SOURCE.nvertices(); i++)
        SUM += current_sim(i);

    SOURCE = tmp;
    return SUM;
}

bool Rigid_cost_function::get_all_neighbours(int index, std::vector<int>& N, int n, const newresampler::Mesh& REF, std::shared_ptr<Neighbourhood>& nbh, MISCMATHS::SpMat<int>& found) {

    bool update = false;

    for (auto j = REF.tIDbegin(n); j != REF.tIDend(n); j++)
    {
        int n0 = REF.get_triangle(*j).get_vertex_no(0),
                n1 = REF.get_triangle(*j).get_vertex_no(1),
                n2 = REF.get_triangle(*j).get_vertex_no(2);

        if((*nbh)(index, 0) != n0
           || (*nbh)(index, 0) != n1
           || (*nbh)(index, 0) != n2)
            update = true;

        if(found.Peek(n0 + 1,1) == 0) { N.push_back(n0); found.Set(n0 + 1,1,1); }
        if(found.Peek(n1 + 1,1) == 0) { N.push_back(n1); found.Set(n1 + 1,1,1); }
        if(found.Peek(n2 + 1,1) == 0) { N.push_back(n2); found.Set(n2 + 1,1,1); }
    }

    return update;
}

newresampler::Mesh Rigid_cost_function::run(){

    if(verbosity) std::cout << "Rigid registration started" << std::endl;

    double Euler1 = 0.0, Euler2 = 0.0, Euler3 = 0.0, RECfinal = 0.0;
    int min_iter = 0, loop = 0;

    double grad_zero = rigid_cost_mesh(Euler1, Euler2, Euler3);
    double mingrad_zero = grad_zero, RECinit = grad_zero;

    while (spacing > 0.05)
    {
        double step = stepsize;
        double per = spacing;

        for (int _iter = 1; _iter <= iters; _iter++)
        {
            //Initially Euler angles are zero
            Euler1 = 0.0; Euler2 = 0.0; Euler3 = 0.0;

            newresampler::Point grad;
            grad.X = (rigid_cost_mesh(Euler1 + per, Euler2, Euler3) - grad_zero) / per;
            grad.Y = (rigid_cost_mesh(Euler1, Euler2 + per, Euler3) - grad_zero) / per;
            grad.Z = (rigid_cost_mesh(Euler1, Euler2, Euler3 + per) - grad_zero) / per;
            grad.normalize();

            Euler1 += step * grad.X;
            Euler2 += step * grad.Y;
            Euler3 += step * grad.Z;

            if (verbosity)
                std::cout << "LOOP " << loop << "\titer " << _iter << "\tper " << per << "\tstep: " << step
                          << "\n grad_zero: " << grad_zero << "\tmingrad_zero: " << mingrad_zero
                          << " (loop*iters)+_iter " << (loop * iters) + _iter << "\tmin_iter " << min_iter << std::endl;


            newresampler::Mesh tmp = SOURCE;

            rotate_in_mesh(Euler1, Euler2, Euler3);
            grad_zero = rigid_cost_mesh(Euler1, Euler2, Euler3);

            if(grad_zero > mingrad_zero)
            {
                mingrad_zero = grad_zero;
                min_iter = (loop * iters) + _iter;
                RECfinal = mingrad_zero;
            }

            if((loop * iters) + _iter - min_iter > 0)
            {
                step *= 0.5;
                SOURCE = tmp;
            }

            if(step < 1e-3) break;
        }
        loop++;
        spacing *= 0.5;
    }

    if(verbosity && (RECfinal != 0.0))
        std::cout << "Rigid registration finished. Affine improvement: " << abs(((RECfinal-RECinit))/RECfinal*100.0) << "%" << std::endl;

    return SOURCE;
}

} //namespace newmeshreg
