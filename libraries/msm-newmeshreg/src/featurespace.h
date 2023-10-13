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
#ifndef NEWMESHREG_FEATURESPACE_H
#define NEWMESHREG_FEATURESPACE_H

#include "newresampler/resampler.h"
#include "reg_tools.h"

#define RAD 100.0
#define EPSILON 1.0E-8

namespace newmeshreg {

class featurespace {

public:
    featurespace(const std::string& datain, const std::string& dataref);
    featurespace(const std::vector<std::string>& datalist);

    //---INITIALISE---//
    newresampler::Mesh initialize(int, std::vector<newresampler::Mesh>&, bool);

    //---SET---//
    inline void set_smoothing_parameters(const std::vector<double>& s) { _sigma_in = s; }
    inline void set_cutthreshold(const std::vector<float>& thr) { _fthreshold = thr; }
    inline void varnorm(bool norm) { _varnorm = norm; }
    inline void is_sparse(bool sp) { _issparse = sp; }
    inline void intensitynormalize(bool norm, bool _exclcut) { _intensitynorm = norm; _cut = _exclcut; }
    inline void set_nthreads(int nthreads) { _nthreads = nthreads; }

    //---GET---//
    inline int get_dim() const { return DATA[0]->Nrows(); }
    inline double get_input_val(int i, int j) const { return DATA[0]->Peek(i, j); }
    inline double get_ref_val(int i, int j) const { return DATA[1]->Peek(i,j); }
    inline std::shared_ptr<newresampler::Mesh> get_input_excl() const { return EXCL[0]; }
    inline std::shared_ptr<newresampler::Mesh> get_reference_excl() const { return EXCL[1]; }
    inline std::shared_ptr<MISCMATHS::BFMatrix> get_input_data() const { return DATA[0]; }
    inline std::shared_ptr<MISCMATHS::BFMatrix> get_reference_data() const { return DATA[1]; }
    inline NEWMAT::Matrix get_data_matrix(int i) const { return DATA[i]->AsMatrix(); }
    inline double get_data_val(int i, int j, int n) const { return DATA[n]->Peek(i,j); }

private:
    std::vector<std::shared_ptr<MISCMATHS::BFMatrix>> DATA; // holds generic BFMATRIX data which can be sparse or full matrices
    std::vector<std::shared_ptr<newresampler::Mesh>> EXCL;  // exclusion masks for binary weighting of the data during resampling

    std::vector<std::string> CMfile_in;  // path to data
    std::vector<double> _sigma_in;  // smoothing parameters for input and reference
    std::vector<float> _fthreshold;

    bool _intensitynorm = false; // will histogram match
    bool _issparse = false;  // notes that data is sparse
    bool _varnorm = false;  // performs online variance normalisation - maybe replace with non online version called during logtransformandnormalise()
    bool _cut = false;
    int _nthreads = 1;

    //---PROCESSING---//
    void variance_normalise(std::shared_ptr<MISCMATHS::BFMatrix>&, std::shared_ptr<newresampler::Mesh>&);
};

} //namespace newmeshreg

#endif //NEWMESHREG_FEATURESPACE_H
