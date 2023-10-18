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
#include "featurespace.h"

namespace newmeshreg {

featurespace::featurespace(const std::string& datain, const std::string& dataref) {
    _sigma_in.resize(2,5.0);
    _fthreshold.resize(2,0.0);
    CMfile_in.push_back(datain);
    CMfile_in.push_back(dataref);
}

featurespace::featurespace(const std::vector<std::string>& datalist) {
    _sigma_in.resize(datalist.size(), 5.0);
    _fthreshold.resize(2,0.0);
    CMfile_in = datalist;
}

newresampler::Mesh featurespace::initialize(int ico, std::vector<newresampler::Mesh>& IN, bool exclude) {

    newresampler::Mesh icotmp;

    if (IN.size() != CMfile_in.size())
        throw MeshregException("featurespace::Initialize do not have the same number of datasets and surface meshes");
    else
        DATA.resize(CMfile_in.size(), std::shared_ptr<MISCMATHS::BFMatrix>());

    if(ico > 0)
    {
        icotmp = newresampler::make_mesh_from_icosa(ico);
        newresampler::recentre(icotmp);
        newresampler::true_rescale(icotmp, RAD);
    }

    for (unsigned int i = 0; i < IN.size(); i++)
    {
        set_data(CMfile_in[i],DATA[i],IN[i],_issparse);

        if (ico == 0) icotmp = IN[i];

        if (exclude || _cut)
            EXCL.push_back(std::make_shared<newresampler::Mesh>(
                    newresampler::create_exclusion(IN[i], _fthreshold[0], _fthreshold[1])));
        else
            EXCL.push_back(std::shared_ptr<newresampler::Mesh>());

        newresampler::Mesh tmp = newresampler::metric_resample(IN[i], icotmp, _nthreads, EXCL[i]);

        if (_sigma_in[i] > 0.0)
            tmp = newresampler::smooth_data(tmp, tmp, _sigma_in[i], _nthreads, EXCL[i]);

        DATA[i] = std::make_shared<MISCMATHS::FullBFMatrix>(tmp.get_pvalues());
    }

    // intensity normalise using histogram matching
    if (_intensitynorm)
        for (unsigned int i = 1; i < IN.size(); i++)
            multivariate_histogram_normalization(*DATA[i], *DATA[0], EXCL[i], EXCL[0], _nthreads);
            // match input data feature distributions to equivalent in ref, rescale all to first feature in reference if _scale is

    if (_varnorm)
        for (unsigned int i = 0; i < IN.size(); i++)
            variance_normalise(DATA[i], EXCL[i]);

    return icotmp;
}

void featurespace::variance_normalise(std::shared_ptr<MISCMATHS::BFMatrix>& DATA, std::shared_ptr<newresampler::Mesh>& EXCL) {

    std::vector<std::vector<double>> _data(DATA->Nrows());

    #pragma omp parallel for num_threads(_nthreads)
    for (unsigned int k = 1; k <= DATA->Nrows(); k++)
        for (unsigned int i = 1; i <= DATA->Ncols(); i++)
            if (!EXCL || EXCL->get_pvalue(i - 1) > 0.0)
                _data[k - 1].push_back(DATA->Peek(k, i));

    std::vector<double> mean(DATA->Nrows(), 0.0),
                        var(DATA->Nrows(), 0.0);

    #pragma omp parallel for num_threads(_nthreads)
    for (unsigned int i = 0; i < _data.size(); ++i)
    {
        for(unsigned int j = 0; j < _data[i].size(); j++)
        {
            double delta = _data[i][j] - mean[i];
            mean[i] += delta / (j + 1);
            var[i] += delta * (_data[i][j] - mean[i]);
        }

        var[i] /= (_data[i].size()-1);

        for(unsigned int j = 0; j < _data[i].size(); ++j)
        {
            _data[i][j] -= mean[i];
            if(var[i] > 0.0) _data[i][j] /= std::sqrt(var[i]);
        }
    }

    #pragma omp parallel for num_threads(_nthreads)
    for(int i = 1; i <= DATA->Nrows(); ++i)
    {
        int data_index = 0;
        for (int j = 1; j <= DATA->Ncols(); ++j)
            if (!EXCL || EXCL->get_pvalue(j - 1) > 0.0)
                DATA->Set(i, j, _data[i - 1][data_index++]);
    }
}

} //namespace newmeshreg
