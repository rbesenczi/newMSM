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
#ifndef NEWMESHREG_SIMILARITIES_H
#define NEWMESHREG_SIMILARITIES_H

#include <miscmaths/bfmatrix.h>

#include "reg_tools.h"

namespace newmeshreg {

class sparsesimkernel{

public:
    sparsesimkernel(): mp(std::make_shared<MISCMATHS::SpMat<double>>()) {}

    //---FOR RIGID---//
    void initialise(int simval);
    inline double peek(unsigned int r, unsigned int c) const { return(mp->Peek(r,c)); }
    void set_simval(int val) { _sim = val; }
    void set_percentile(double val) { percentile = val; }
    void set_input(std::shared_ptr<MISCMATHS::BFMatrix> in ){ m_A = in; }
    void set_reference(std::shared_ptr<MISCMATHS::BFMatrix> ref){ m_B = ref; }
    void set_neighbourhood(std::shared_ptr<Neighbourhood>& n) { nbh = n; }
    void resize(unsigned int m, unsigned int n) { mp = std::make_shared<MISCMATHS::SpMat<double>>(m, n); }
    void calculate_sim_column_nbh(int);

    //---FOR DISCRETE---//
    double get_sim_for_min(const std::vector<double>& input, const std::vector<double>& reference, const std::vector<double>& weights = std::vector<double>()) {
        if(_sim == 1)
            return SSD(input,reference,weights);
        if(_sim == 2)
            return 1 - ( 1 + corr(input,reference,weights)) * 0.5;
        if(_sim == 4)
            return DICE(input,reference);
        if(_sim == 5)
            return genDICE(input,reference);
        throw MeshregException("Unknown similarity metric");
    }

private:
    // DATA
    std::shared_ptr<MISCMATHS::BFMatrix> m_A, m_B;
    std::shared_ptr<MISCMATHS::SpMat<double>> mp;
    std::shared_ptr<Neighbourhood> nbh;
    NEWMAT::RowVector _rmeanA, _rmeanB;

    int _sim = 1;
    double percentile = 0.75;

    //---FOR RIGID---//
    double corr(int, int);
    double SSD(int, int);
    NEWMAT::RowVector meanvector(const MISCMATHS::BFMatrix &); // for correlation measure

    //---FOR DISCRETE---//
    double corr(const std::vector<double>& A, const std::vector<double>& B, const std::vector<double>& weights);
    double corr(const std::vector<double>& A, const std::vector<double>& B);
    double SSD(const std::vector<double>& A, const std::vector<double>& B, const std::vector<double>& weights);
    double SSD(const std::vector<double>& A, const std::vector<double>& B);
    double DICE(const std::vector<double>& A, const std::vector<double>& B);
    double genDICE(const std::vector<double>& A, const std::vector<double>& B);
};

} //namespace newmeshreg

#endif //NEWMESHREG_SIMILARITIES_H
