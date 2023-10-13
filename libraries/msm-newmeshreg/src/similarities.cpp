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
#include "similarities.h"

namespace newmeshreg {

//--------------- FOR RIGID COST FUNCTIONS ---------------//
void sparsesimkernel::Initialize(int simval) {  // for rigid

    _sim = simval;

    if (m_A == nullptr) throw newmeshreg::MeshregException("SIMILARITIES:: Connectivity matrices have not been initliased");

    _rmeanA = meanvector(*m_A);
    if(m_A->Nrows() == 1) _meanA = _rmeanA(1);
    if(m_B != nullptr)
        _rmeanB = meanvector(*m_B);
    else
        _rmeanB = _rmeanA;

    if(m_B->Nrows() == 1) _meanB = _rmeanB(1);
}

void sparsesimkernel::calculate_sim_column_nbh(int ind) {

    for (int j = 0; j < nbh->nrows(ind); j++)
        if ((*nbh)(ind, j) != 0)
            switch (_sim)
            {
                case 1:
                    mp->Set((*nbh)(ind, j) + 1, ind + 1, -SSD(ind + 1, (*nbh)(ind, j) + 1));
                    break;
                case 2:
                    mp->Set((*nbh)(ind, j) + 1, ind + 1, corr(ind + 1, (*nbh)(ind, j) + 1));
                    break;
            }
}

double sparsesimkernel::corr(int i, int j) {

    int num = 0, numB = 0;
    double prod = 0.0, varA = 0.0, varB = 0.0;

    std::shared_ptr<MISCMATHS::SparseBFMatrix<double>> ptr = std::dynamic_pointer_cast<MISCMATHS::SparseBFMatrix<double> >(m_A);

    double Bzerooffset = (0.0-_rmeanB(j));  // result for all zero values of sparse mat
    double Azerooffset = (0.0-_rmeanA(i));

    for (MISCMATHS::BFMatrixColumnIterator it = m_A->begin(i); it != m_A->end(i); it++)
    {
        prod += (*it-_rmeanA(i))*(m_B->Peek(it.Row(),j)-_rmeanB(j));
        varA += (*it-_rmeanA(i))*(*it-_rmeanA(i));
        varB += (m_B->Peek(it.Row(),j)-_rmeanB(j))*(m_B->Peek(it.Row(),j)-_rmeanB(j));
        num++;
    }

    if(ptr)
    { //if sparse run all for all rows where A had zero values

        varA += (m_A->Nrows()-num)*Azerooffset*Azerooffset; // for rows where A has no values
        for (MISCMATHS::BFMatrixColumnIterator it = m_B->begin(j); it != m_B->end(j); it++)
        {
            if(m_A->Peek(it.Row(),i) == 0)
            {
                prod +=Azerooffset*(*it-_rmeanB(j));
                num++;
                varB+=(*it-_rmeanB(j))*(*it-_rmeanB(j));
            }
            numB++;
        }
        varB += (m_A->Nrows()-numB)*Bzerooffset*Bzerooffset;    // for rows where B has no values
        prod += (2*m_A->Nrows()-num)*Azerooffset*Bzerooffset;  // for rows where A&B have no values
    }

    if (varA == 0.0 || varB == 0.0) return 0.0;
    else return prod/(sqrt(varA)*sqrt(varB));
}

double sparsesimkernel::SSD(int i, int j) {

    double prod = 0.0;
    if(m_B == nullptr)  m_B = m_A;

    std::shared_ptr<MISCMATHS::SparseBFMatrix<double>> ptr = std::dynamic_pointer_cast<MISCMATHS::SparseBFMatrix<double>>(m_A);

    for (MISCMATHS::BFMatrixColumnIterator it = m_A->begin(i); it != m_A->end(i); it++)
        prod += (*it-m_B->Peek(it.Row(),j))*(*it-m_B->Peek(it.Row(),j));

    if (ptr)
        for (MISCMATHS::BFMatrixColumnIterator it = m_B->begin(j); it != m_B->end(j); it++)
            if (m_A->Peek(it.Row(), i) == 0)
                prod += (*it) * (*it);

    return sqrt(prod)/m_A->Nrows();
}

NEWMAT::RowVector sparsesimkernel::meanvector(const MISCMATHS::BFMatrix& fdt_matrix) {

    NEWMAT::RowVector mean(fdt_matrix.Ncols());
    mean = 0;

    if (fdt_matrix.Nrows() == 1)
    {
        double sum = 0;
        for (unsigned int i = 0; i < fdt_matrix.Ncols(); i++)
            sum = sum + fdt_matrix.Peek(1, i + 1);
        for (unsigned int i = 0; i < fdt_matrix.Ncols(); i++)
            mean(i + 1) = sum / fdt_matrix.Ncols();
    }
    else
        for (unsigned int i = 0; i < fdt_matrix.Ncols(); i++)
        {
            double sum = 0;
            for (unsigned int j = 0; j < fdt_matrix.Nrows(); j++)
                sum = sum + fdt_matrix.Peek(j + 1, i + 1);
            mean(i + 1) = sum / fdt_matrix.Nrows();
        }

    return mean;
}

//--------------- FOR DISCRETE COST FUNCTIONS ---------------//
// for the case where we are working on vector data and we have no knowledge of the full data m_A (currently used in discrete opt)
void sparsesimkernel::initialize(const std::vector<double>& inputdata, const std::vector<double>& refdata, const std::vector<double>& weights) {

    _meanA = 0.0;
    _meanB = 0.0;
    double sum = 0.0;

    for(unsigned int i = 0; i < inputdata.size(); i++)
    {
        double varwght = 1.0;
        if(weights.size() == inputdata.size()) varwght = weights[i];

        _meanA += varwght * inputdata[i];
        _meanB += varwght * refdata[i];
        sum += varwght;
    }

    if(sum > 0.0)
    {
        _meanA/=sum;
        _meanB/=sum;
    }
}

double sparsesimkernel::corr(const std::vector<double>& A, const std::vector<double>& B, const std::vector<double>& weights) {

    double prod = 0.0, varA = 0.0, varB = 0.0;
    double sum = 0.0;

    initialize(A, B, weights);

    for (unsigned int s = 0; s < A.size(); s++)
    {
        double varwght = 1.0;

        if (weights.size() == A.size()) varwght = weights[s];

        prod += varwght * (A[s] - _meanA) * (B[s] - _meanB);
        varA += varwght * (A[s] - _meanA) * (A[s] - _meanA);
        varB += varwght * (B[s] - _meanB) * (B[s] - _meanB);
        sum += varwght;
    }

    if(sum > 0)
    {
      prod /= sum;
      varA /= sum;
      varB /= sum;
    }

    if (varA == 0.0 || varB == 0.0) return 0.0;
    else return prod / (sqrt(varA) * sqrt(varB));
}

double sparsesimkernel::corr(const std::map<int,float> &A, const std::map<int,float> &B) {

    double prod = 0.0, varA = 0.0, varB = 0.0, MAPA_MEAN = 0.0, MAPB_MEAN = 0.0;
    int ind = 0;

    // estimate mean
    for (const auto &iter: A)
        if (B.find(iter.first) != B.end()) {
            MAPA_MEAN += iter.second;
            MAPB_MEAN += B.find(iter.first)->second;
            ind++;
        }

    if (ind == 0) return 1e6;
    else {
        MAPA_MEAN /= ind;
        MAPB_MEAN /= ind;

        for (const auto &iter: A) {
            if (B.find(iter.first) != B.end()) {
                prod += (iter.second - MAPA_MEAN) * (B.find(iter.first)->second - MAPB_MEAN);
                varA += (iter.second - MAPA_MEAN) * (iter.second - MAPA_MEAN);
                varB += (B.find(iter.first)->second - MAPB_MEAN) * (B.find(iter.first)->second - MAPB_MEAN);
            }
        }

        prod /= ind;
        varA /= ind;
        varB /= ind;

        if (varA == 0.0 || varB == 0.0) return 0.0;
        else return (1 - (1 + (prod / (sqrt(varA) * sqrt(varB)))) * 0.5);
    }
}
} //namespace newmeshreg
