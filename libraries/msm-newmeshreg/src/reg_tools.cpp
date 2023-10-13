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
#include "reg_tools.h"

namespace newmeshreg {

const char* MeshregException::what() const noexcept {
    std::cout << errmesg << std::endl;
    return errmesg;
}

void Neighbourhood::update(const newresampler::Mesh& source, const newresampler::Mesh& target, double ang, int numthreads) {

    angsep = ang;
    neighbours.clear();
    neighbours.resize(source.nvertices());

    #pragma omp parallel for num_threads(numthreads)
    for(int index = 0; index < source.nvertices(); ++index)
    {
        newresampler::Point cr = source.get_coord(index);
        cr.normalize();
        std::vector<std::pair<double, int>> cr_neighbours;
        for(int n = 0; n < target.nvertices(); ++n)
        {
            newresampler::Point actual = target.get_coord(n);
            actual.normalize();
            if ((actual | cr) >= cos(angsep))
                cr_neighbours.emplace_back(std::make_pair((actual - cr).norm(), n));
        }

        std::sort(cr_neighbours.begin(), cr_neighbours.end(),
                  [](const auto& lhs, const auto& rhs) -> bool { return lhs.first < rhs.first; });

        for(const auto& n : cr_neighbours)
            neighbours[index].emplace_back(n.second);
    }
}

bool get_all_neighbours(int index, std::vector<int>& N, const newresampler::Point& point, int n,
                        const newresampler::Mesh& REF, std::shared_ptr<Neighbourhood>& nbh, MISCMATHS::SpMat<int>& found) {

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

void computeNormal2EdgeOfTriangle(const newresampler::Point& v0, const newresampler::Point& v1, const newresampler::Point& v2, newresampler::Point& norm2edge) {

    newresampler::Point s1 = v2 - v0, s2 = v1 - v0;

    if(s1.norm() > 1e-10)  s1.normalize(); else s1 = newresampler::Point(0,0,0);
    if(s2.norm() > 1e-10)  s2.normalize(); else s2 = newresampler::Point(0,0,0);

    newresampler::Point norm2triangle = s1 * s2;

    if(norm2triangle.norm() > 1e-10)
        norm2triangle.normalize();
    else
        norm2triangle = newresampler::Point(0,0,0);

    norm2edge = s2 * norm2triangle;
    //Then take cross product of norm and the edge v1-v0 to get vector perpendicular to v1-v0

    if ((s1 | norm2edge) < 0) /* first edge is in wrong direction, flip it. */
        norm2edge = norm2edge * -1;
}

newresampler::Point computeGradientOfBarycentricTriangle(const newresampler::Point& v0, const newresampler::Point& v1, const newresampler::Point& v2) {
// Find the gradient of triangle area associated with edge v1v0 (end moving point v2) points away from v1v0 with magnitude equal to half the length of v1v0
// used for unfolding

    newresampler::Point v0v1, norm2edge;
    double base;

    computeNormal2EdgeOfTriangle(v0, v1, v2, norm2edge);

    v0v1 = v1 - v0;
    base = v0v1.norm();

    return norm2edge * 0.5 * base;
}

newresampler::Point spatialgradient(int index, const newresampler::Mesh& SOURCE) {
// spatial gradient for vertex triangular mesh
    newresampler::Point grad, ci = SOURCE.get_coord(index);

    for (auto j = SOURCE.tIDbegin(index); j != SOURCE.tIDend(index); j++)
    {
        newresampler::Point v0 = SOURCE.get_triangle_vertex(*j, 0),
                v1 = SOURCE.get_triangle_vertex(*j, 1),
                v2 = SOURCE.get_triangle_vertex(*j, 2),
                dA;

        if ((ci - v0).norm() == 0)
            dA = computeGradientOfBarycentricTriangle(v1, v2, v0);
        else if ((ci - v1).norm() == 0)
            dA = computeGradientOfBarycentricTriangle(v2, v0, v1);
        else
            dA = computeGradientOfBarycentricTriangle(v0, v1, v2);

        grad = grad + dA;
    }
    return grad;
}

bool check_for_intersections(int ind, newresampler::Mesh& IN) {

    newresampler::Point N = IN.get_triangle_from_vertex(ind, 0).normal();

    for (auto j = IN.tIDbegin(ind); j != IN.tIDend(ind); j++)
    {
        newresampler::Point N2 = IN.get_triangle(*j).normal();
        if ((N | N2) <= 0.5) return true;
    }

    return false;
}

void unfold(newresampler::Mesh& SOURCE, bool verbosity) {

    std::vector<int> foldedvertices;
    std::vector<newresampler::Point> foldinggradients;
    bool folded = true;
    int it = 0;
    const double RAD = 100;

    while (folded)
    {
        foldedvertices.clear();
        foldinggradients.clear();

        for (int i = 0; i < SOURCE.nvertices(); i++)
            if (check_for_intersections(i, SOURCE))
                foldedvertices.push_back(i);

        if (foldedvertices.empty()) return;

        else if (it % 100 == 0 && verbosity)
            std::cout << "Mesh is folded, total folded vertices: " << foldedvertices.size() << " iter: " << it << std::endl;

        for (int foldedvertice: foldedvertices)
            foldinggradients.emplace_back(spatialgradient(foldedvertice, SOURCE));
            // estimates gradient of triangle area for all vertices that are detected as folded

        for (unsigned int i = 0; i < foldedvertices.size(); i++)
        {
            double current_stepsize = 1.0;
            newresampler::Point ci = SOURCE.get_coord(foldedvertices[i]);
            newresampler::Point grad = foldinggradients[i];
            newresampler::Point pp;

            do
            {
                pp = ci - grad * current_stepsize; // gradient points in direction of increasing area therefore opposite to desired direction
                pp.normalize();
                SOURCE.set_coord(foldedvertices[i],pp * RAD);
                current_stepsize *= 0.5;
            }
            while(check_for_intersections(foldedvertices[i], SOURCE) && current_stepsize > 1e-3);

            SOURCE.set_coord(foldedvertices[i],pp * RAD);
        }
        it++;
        if(it == 1000) break;
    }
}

NEWMAT::Matrix get_coordinate_transformation(double dNdT1,double dNdT2, NEWMAT::ColumnVector& Norm) {

    newresampler::Point G1(1, 0, dNdT1);
    newresampler::Point G2(0, 1, dNdT2);
    newresampler::Point G3 = G1 * G2;

    G3 = G3 / sqrt((G3|G3));

    NEWMAT::Matrix G(3,3);

    G(1,1) = G1.X;
    G(1,2) = G2.X;
    G(1,3) = G3.X;
    G(2,1) = G1.Y;
    G(2,2) = G2.Y;
    G(2,3) = G3.Y;
    G(3,1) = G1.Z;
    G(3,2) = G2.Z;
    G(3,3) = G3.Z;
    Norm(1) = G3.X;
    Norm(2) = G3.Y;
    Norm(3) = G3.Z;

    return (G.i()).t();
}

newresampler::Tangs calculate_tangs(int ind, const newresampler::Mesh& SPH_in) {

    newresampler::Tangs T;
    double mag;
    newresampler::Point a = SPH_in.local_normal(ind);
    //newresampler::Point tmp = SPH_in.get_coord(ind);

    if((a|SPH_in.get_coord(ind)) < 0) a = a * -1;

    if(abs(a.X) >= abs(a.Y) && abs(a.X) >= abs(a.Z))
    {
        mag = sqrt(a.Z*a.Z + a.Y*a.Y);
        if(mag == 0)
        {
            T.e1.X = 0;
            T.e1.Y = 0;
            T.e1.Z = 1;
        }
        else
        {
            T.e1.X = 0;
            T.e1.Y = -a.Z /mag;
            T.e1.Z = a.Y /mag;
        }
    }
    else if(abs(a.Y) >= abs(a.X) && abs(a.Y) >= abs(a.Z))
    {
        mag = sqrt(a.Z*a.Z + a.X*a.X);
        if(mag == 0)
        {
            T.e1.X = 0;
            T.e1.Y = 0;
            T.e1.Z = 1;
        }
        else
        {
            T.e1.X = -a.Z/mag;
            T.e1.Y = 0;
            T.e1.Z = a.X /mag;
        }
    }
    else
    {
        mag = sqrt(a.Y*a.Y + a.X*a.X);
        if(mag == 0)
        {
            T.e1.X = 1;
            T.e1.Y = 0;
            T.e1.Z = 0;
        }
        else
        {
            T.e1.X = -a.Y/mag;
            T.e1.Y = a.X/mag;
            T.e1.Z = 0;
        }
    }
    T.e2 = a * T.e1;
    T.e2.normalize();

    return T;
}

newresampler::Tangs calculate_tri(const newresampler::Point& a) {

    newresampler::Tangs T;
    newresampler::Point b(1.0f,0.0f,0.0f);
    newresampler::Point c = a * b;

    double len = c.X * c.X + c.Y * c.Y + c.Z * c.Z;

    if (len == 0.0)
    {
        /* the vector b was parallel to a */
        b.X = 0.0f;
        b.Y = 1.0f;
        b.Z = 0.0f;
        c = a * b;
        len = c.X * c.X + c.Y * c.Y + c.Z * c.Z;
    }

    /* normalize */
    len = std::sqrt(len);

    if (len == 0.0)
    {
        std::cout << "Tangs Tangent:: first tangent vector of length zero at vertex: " << std::endl;
        len = 1;
    }

    T.e1.X = c.X / len ;
    T.e1.Y = c.Y / len ;
    T.e1.Z = c.Z / len ;

    b = a * c;
    /* normalize */
    len = std::sqrt(b.X * b.X + b.Y * b.Y + b.Z * b.Z);

    if (len == 0)
    {
        std::cout << "Tangs Tangent::second tangent vector of length zero at vertex: " << std::endl;
        len = 1 ;
    }

    T.e2.X = b.X / len ;
    T.e2.Y = b.Y / len ;
    T.e2.Z = b.Z / len ;

    return T;
}

newresampler::Tangs calculate_tri(int ind, const newresampler::Mesh& SPH_in) {

    newresampler::Tangs T;
    newresampler::Point a = SPH_in.get_triangle_normal(ind);
    newresampler::Point b(1.0f,0.0f,0.0f);
    newresampler::Point c = a * b;

    double len = c.X * c.X + c.Y * c.Y + c.Z * c.Z;

    if (len == 0.0)
    {
        /* the vector b was parallel to a */
        b.X = 0.0f;
        b.Y = 1.0f;
        b.Z = 0.0f;
        c = a * b;
        len = c.X * c.X + c.Y * c.Y + c.Z * c.Z;
    }

    /* normalize */
    len = std::sqrt(len);

    if (len == 0.0)
    {
        std::cout << "Tangs Tangent:: first tangent vector of length zero at vertex: " << ind << std::endl;
        len = 1;
    }

    T.e1.X = c.X / len ;
    T.e1.Y = c.Y / len ;
    T.e1.Z = c.Z / len ;

    b = a * c;

    /* normalize */
    len = std::sqrt(b.X * b.X + b.Y * b.Y + b.Z * b.Z);

    if (len == 0)
    {
        std::cout << "Tangs Tangent::second tangent vector of length zero at vertex: " << ind << std::endl;
        len = 1;
    }

    T.e2.X = b.X / len ;
    T.e2.Y = b.Y / len ;
    T.e2.Z = b.Z / len ;

    return T;
}

NEWMAT::ColumnVector calculate_strains(int index, const std::vector<int>& kept, const newresampler::Mesh& orig, const newresampler::Mesh& final, const std::shared_ptr<NEWMAT::Matrix>& PrincipalStretches){

    NEWMAT::ColumnVector STRAINS(4); STRAINS = 0;
    newresampler::Point Normal_O = orig.get_normal(index);
    NEWMAT::Matrix TRANS(3,3), a, b, c, d;
    //double T1_coord = 0.0, T2_coord = 0.0;
    int maxind, minind;
    newresampler::Tangs T = calculate_tangs(index, orig);
    TRANS = form_matrix_from_points(T.e1, T.e2, Normal_O, true);
    NEWMAT::Matrix pseudo_inv, pseudo_V, pseudo_U;
    NEWMAT::DiagonalMatrix pseudo_D;
    NEWMAT::Matrix alpha(kept.size(),6),
                   N(kept.size(),1),
                   n(kept.size(),1),
                   t1(kept.size(),1),
                   t2(kept.size(),1);
    alpha = 0; N = 0; n = 0; t1 = 0; t2 = 0;

    for(int j = 0; j < (int)kept.size(); j++)
    {
        double T1_coord = 0.0, T2_coord = 0.0;
        newresampler::Point tmp = orig.get_coord(kept[j]) - orig.get_coord(index);
        project_point(tmp, T, T1_coord, T2_coord);
        N(j+1,1) = tmp | Normal_O;

        // fit poly
        alpha(j+1,2) = T1_coord;
        alpha(j+1,3) = T2_coord;
        alpha(j+1,4) = 0.5 * T1_coord * T1_coord;
        alpha(j+1,5) = 0.5 * T2_coord * T2_coord;
        alpha(j+1,6) = T1_coord * T2_coord;

        // get transformed coords
        tmp = final.get_coord(kept[j]) - final.get_coord(index);
        project_point(tmp, T, T1_coord, T2_coord);

        n(j + 1,1) = tmp | Normal_O;
        t1(j + 1,1) = T1_coord;
        t2(j + 1,1) = T2_coord;
    }

    SVD(alpha, pseudo_D, pseudo_U, pseudo_V);

    for (int i = 1; i <= pseudo_D.Nrows(); i++)
        if (pseudo_D(i) != 0)
            pseudo_D(i) = 1 / pseudo_D(i);

    pseudo_inv = pseudo_V * pseudo_D * pseudo_U.t();

    // get coeffieicnts for different fits
    a = pseudo_inv * N;
    b = pseudo_inv * t1;
    c = pseudo_inv * t2;
    d = pseudo_inv * n;

    double dNdT1, dNdT2, dt1dT1, dt1dT2, dt2dT1, dt2dT2, dndT1, dndT2;
    dNdT1 = a(2,1); dNdT2 = a(3,1);
    dt1dT1 = b(2,1); dt1dT2 = b(3,1);
    dt2dT1 = c(2,1); dt2dT2 = c(3,1);
    dndT1 = d(2,1); dndT2 = d(3,1);

    NEWMAT::ColumnVector G3(3);
    NEWMAT::Matrix G_cont = get_coordinate_transformation(dNdT1, dNdT2, G3);

    newresampler::Point g1(dt1dT1,dt2dT1,dndT1);
    newresampler::Point g2(dt1dT2,dt2dT2,dndT2);
    newresampler::Point g3 = g1 * g2; g3 = g3 / sqrt((g3|g3));

    NEWMAT::Matrix g(3,3);
    g(1,1)=g1.X; g(1,2)=g2.X; g(1,3)=g3.X;
    g(2,1)=g1.Y; g(2,2)=g2.Y; g(2,3)=g3.Y;
    g(3,1)=g1.Z; g(3,2)=g2.Z; g(3,3)=g3.Z;

    NEWMAT::Matrix F = g * G_cont.t();
    NEWMAT::Matrix C = F.t() * F;

    NEWMAT::DiagonalMatrix Omega;
    NEWMAT::Matrix U, Umax, Umin;
    SVD(C,Omega,U);

    NEWMAT::Matrix mm = G3.t() * U;

    if (abs(mm(1, 1)) >= abs(mm(1, 2)) && abs(mm(1, 1)) >= abs(mm(1, 3)))
        if (sqrt(Omega(2)) > sqrt(Omega(3)))
        {
            maxind = 2;
            minind = 3;
        }
        else
        {
            maxind = 3;
            minind = 2;
        }

    else if (abs(mm(1, 2)) >= abs(mm(1, 1)) && abs(mm(1, 2)) >= abs(mm(1, 3)))
        if (sqrt(Omega(1)) > sqrt(Omega(3)))
        {
            maxind = 1;
            minind = 3;
        }
        else
        {
            maxind = 3;
            minind = 1;
        }
    else
    {
        if (sqrt(Omega(1)) > sqrt(Omega(2)))
        {
            maxind = 1;
            minind = 2;
        }
        else
        {
            maxind = 2;
            minind = 1;
        }
    }

    STRAINS(1) = sqrt(Omega(maxind));
    STRAINS(2) = sqrt(Omega(minind));
    Umax = TRANS.i() * U.SubMatrix(1,3,maxind,maxind);
    Umin = TRANS.i() * U.SubMatrix(1,3,minind,minind);

    if(PrincipalStretches)
    {
        (*PrincipalStretches)(index+1,1) = Umax(1,1);
        (*PrincipalStretches)(index+1,2) = Umax(2,1);
        (*PrincipalStretches)(index+1,3) = Umax(3,1);
        (*PrincipalStretches)(index+1,4) = Umin(1,1);
        (*PrincipalStretches)(index+1,5) = Umin(2,1);
        (*PrincipalStretches)(index+1,6) = Umin(3,1);
    }

    STRAINS(3) = 0.5 * (STRAINS(1)*STRAINS(1) - 1);
    STRAINS(4) = 0.5 * (STRAINS(2)*STRAINS(2) - 1);

    return STRAINS;
}

newresampler::Mesh calculate_strains(double fit_radius, const newresampler::Mesh& orig, const newresampler::Mesh &final,
                                     int numthreads, const std::shared_ptr<NEWMAT::Matrix>& PrincipalStretches){

    newresampler::Mesh strain = final;
    NEWMAT::Matrix strains(4, orig.nvertices());

    if(PrincipalStretches) PrincipalStretches->ReSize(orig.nvertices(), 6);

    #pragma omp parallel for num_threads(numthreads)
    for (int index = 0; index < orig.nvertices(); index++)
    {
        std::vector<int> kept;

        double fit_temp = fit_radius;
        while (kept.size() <= 8)
        {
            kept.clear();
            for(int j = 0; j < orig.nvertices(); j++)
            {
                double dist = (orig.get_coord(index) - orig.get_coord(j)).norm();
                // reject points with normals in opposite directions
                double dir_chk1 = orig.get_normal(j) | orig.get_normal(index);
                //double dir_chk2 = 1;

                if(dist <= fit_temp && dir_chk1 >= 0/* && dir_chk2 >= 0*/)
                    kept.push_back(j);
            }

            if(kept.size() > 8)
            {
                NEWMAT::ColumnVector node_strain = calculate_strains(index, kept, orig, final, PrincipalStretches);
                strains(1, index + 1) = node_strain(1);
                strains(2, index + 1) = node_strain(2);
                strains(3, index + 1) = node_strain(3);
                strains(4, index + 1) = node_strain(4);
            }
            else
                fit_temp += 0.5;
        }
    }

    strain.set_pvalues(strains);

    return strain;
}

double triangle_strain(const NEWMAT::Matrix& AA, const NEWMAT::Matrix & BB, double MU, double KAPPA, const std::shared_ptr<NEWMAT::ColumnVector>& strains, double k_exp) {

    double c0,c1,c2,c3,c4,c5,c0c,c1c,c2c,c3c,c4c,c5c,I1,I3,J,W;
    NEWMAT::Matrix Edges(2,2), edges(2,2);
    NEWMAT::Matrix F, F3D(3,3), F3D_2;

    F3D(3,1) = 0; F3D(3,2) = 0; F3D(3,3) = 1;

    // all t0 distances
    c0=AA(2,1)-AA(1,1); //x2-x1 t0
    c1=AA(2,2)-AA(1,2); //y2-y1
    c2=AA(3,1)-AA(2,1); //x3-x2
    c3=AA(3,2)-AA(2,2); // y3-y2
    c4=AA(3,1)-AA(1,1); // x3-x1
    c5=AA(3,2)-AA(1,2); // y3-y1

    //all tf distances
    c0c=BB(2,1)-BB(1,1); //x2-x1 tf
    c1c=BB(2,2)-BB(1,2); // y2-y1
    c2c=BB(3,1)-BB(2,1); // x3-x2
    c3c=BB(3,2)-BB(2,2); // y3-y2
    c4c=BB(3,1)-BB(1,1);  // x3-x1
    c5c=BB(3,2)-BB(1,2); // y3-y1

    Edges(1,1)=c0; Edges(1,2)=c4; Edges(2,1)=c1; Edges(2,2)=c5;
    edges(1,1)=c0c; edges(1,2)=c4c; edges(2,1)=c1c; edges(2,2)=c5c;

    F = edges * Edges.i();

    F3D(1,1)=F(1,1); F3D(1,2)=F(1,2); F3D(1,3)=0;
    F3D(2,1)=F(2,1); F3D(2,2)=F(2,2); F3D(2,3)=0;

    F3D_2 = F3D.t() * F3D;
    I1 = F3D_2.Trace();
    I3 = F3D_2.Determinant();

    J = sqrt(I3);
    double I1st_new = (I1 - 1.0) / J;
    double R;
    if (I1st_new <= 2)
        R = 1.0;
    else
        R = 0.5 * (I1st_new + sqrt(I1st_new * I1st_new - 4));//convert I1* to (major strain) / (minor strain)

    double Rshared = pow(R, k_exp), Jshared = pow(J, k_exp);//could use different exponents, but would be quite strange
    W = 0.5 * (MU * (Rshared + 1.0 / Rshared - 2) + KAPPA * (Jshared + 1.0 / Jshared - 2));

    // calculating prinicipal strains (for testing)
    if(strains)
    {
        double a11,a21,a31,a12,a22,a32,a13,a23,a33;
        double A, A11,A21,A31,A12,A22,A32,A13,A23,A33;
        double B1,B2,B3;
        double e11,e22,e12,X,Y;
        //2*dx^2, then 2*dy^2, then 2*dx*dy
        a11=2.*c0*c0 ;
        a21=2.*c2*c2 ;
        a31=2.*c4*c4;
        a12=2.*c1*c1 ;
        a22=2.*c3*c3 ;
        a32=2.*c5*c5;
        a13=4.*c0*c1 ;
        a23=4.*c2*c3;
        a33=4.*c4*c5;

        A=a11*a22*a33+a12*a23*a31+a13*a21*a32-a13*a22*a31-a23*a32*a11-a33*a21*a12;

        A11=(a22*a33-a32*a23)/A;
        A12=(a13*a32-a12*a33)/A;
        A13=(a12*a23-a13*a22)/A;
        A21=(a23*a31-a21*a33)/A;
        A22=(a11*a33-a13*a31)/A;
        A23=(a13*a21-a11*a23)/A;
        A31=(a21*a32-a22*a31)/A;
        A32=(a12*a31-a11*a32)/A;
        A33=(a11*a22-a12*a21)/A;

        // deformed distances between points
        B1=c0c*c0c+c1c*c1c-c0*c0-c1*c1;
        B2=c2c*c2c+c3c*c3c-c2*c2-c3*c3;
        B3=c4c*c4c+c5c*c5c-c4*c4-c5*c5;

        // strains wrt x,y coords
        e11=B1*A11+B2*A12+B3*A13;
        e22=B1*A21+B2*A22+B3*A23;
        e12=B1*A31+B2*A32+B3*A33;

        X=e11+e22 ;
        Y=e11-e22;

        (*strains)(1)=X/2+sqrt((Y/2)*(Y/2)+(e12)*(e12));
        (*strains)(2)=X/2-sqrt((Y/2)*(Y/2)+(e12)*(e12));
    }

    return W;
}

double calculate_triangular_strain(int index, const newresampler::Mesh& ORIG, const newresampler::Mesh& FINAL, double mu, double kappa,
                                   const std::shared_ptr<NEWMAT::ColumnVector>& indexSTRAINS, double k_exp) {

    NEWMAT::Matrix ORIG3D(3,3), FINAL3D(3,3),
                   ORIG2D(3,3), FINAL2D(3,3),
                   TRANS, TRANS2;
    newresampler::Point Normal_O = ORIG.get_triangle_normal(index),
          Normal_F = FINAL.get_triangle_normal(index);

    newresampler::Tangs T = calculate_tri(index,ORIG);
    newresampler::Tangs T_trans = calculate_tri(index,FINAL);
    double W = 0.0;
    NEWMAT::Matrix TMP;
    NEWMAT::DiagonalMatrix D, D2;

    TRANS = form_matrix_from_points(T.e1, T.e2, Normal_O);
    if(TRANS.Determinant() < 0)
    {
        TMP = TRANS;
        TMP(1,1) = TRANS(1,2); TMP(2,1) = TRANS(2,2); TMP(3,1) = TRANS(3,2);
        TMP(1,2) = TRANS(1,1); TMP(2,2) = TRANS(2,1); TMP(3,2) = TRANS(3,1);
        TRANS = TMP;
    }

    TRANS2 = form_matrix_from_points(T_trans.e1,T_trans.e2,Normal_F);
    if(TRANS.Determinant() < 0)
    {
        TMP = TRANS2;
        TMP(1,1) = TRANS2(1,2); TMP(2,1) = TRANS2(2,2); TMP(3,1) = TRANS2(3,2);
        TMP(1,2) = TRANS2(1,1); TMP(2,2) = TRANS2(2,1); TMP(3,2) = TRANS2(3,1);
        TRANS2 = TMP;
    }

    for(int i = 0; i < 3; i++)
    {
        newresampler::Point vertex = ORIG.get_triangle_vertex(index, i);
        ORIG3D(i+1,1) = vertex.X; ORIG3D(i+1,2) = vertex.Y; ORIG3D(i+1,3) = vertex.Z;

        vertex = FINAL.get_triangle_vertex(index, i);
        FINAL3D(i+1,1) = vertex.X;  FINAL3D(i+1,2) = vertex.Y;  FINAL3D(i+1,3) = vertex.Z;
    }

    ORIG2D = ORIG3D * TRANS;
    FINAL2D = FINAL3D * TRANS2;

    W = triangle_strain(ORIG2D, FINAL2D, mu, kappa, indexSTRAINS, k_exp);

    return W;
}

double calculate_triangular_strain(const newresampler::Triangle& ORIG_tr, const newresampler::Triangle& FINAL_tr, double mu, double kappa,
                                   const std::shared_ptr<NEWMAT::ColumnVector>& indexSTRAINS, double k_exp) {

    newresampler::Point Normal_O = ORIG_tr.normal(), Normal_F = FINAL_tr.normal();
    NEWMAT::Matrix ORIG3D(3,3), FINAL3D(3,3),
                   ORIG2D(3,3), FINAL2D(3,3),
                   TRANS, TRANS2;

    newresampler::Tangs T = calculate_tri(Normal_O);
    newresampler::Tangs T_trans = calculate_tri(Normal_F);
    double W;
    NEWMAT::Matrix TMP;
    NEWMAT::DiagonalMatrix D, D2;

    TRANS = form_matrix_from_points(T.e1,T.e2,Normal_O);
    if(TRANS.Determinant() < 0)
    {
        TMP = TRANS;
        TMP(1,1)=TRANS(1,2); TMP(2,1)=TRANS(2,2); TMP(3,1)=TRANS(3,2);
        TMP(1,2)=TRANS(1,1); TMP(2,2)=TRANS(2,1); TMP(3,2)=TRANS(3,1);
        TRANS = TMP;
    }
    TRANS2 = form_matrix_from_points(T_trans.e1,T_trans.e2,Normal_F);
    if(TRANS.Determinant() < 0)
    {
        TMP = TRANS2;
        TMP(1,1)=TRANS2(1,2); TMP(2,1)=TRANS2(2,2); TMP(3,1)=TRANS2(3,2);
        TMP(1,2)=TRANS2(1,1); TMP(2,2)=TRANS2(2,1); TMP(3,2)=TRANS2(3,1);
        TRANS2=TMP;
    }

    for(int i=0;i<3;i++){
        newresampler::Point vertex = ORIG_tr.get_vertex_coord(i);
        ORIG3D(i+1,1) = vertex.X; ORIG3D(i+1,2) = vertex.Y; ORIG3D(i+1,3) = vertex.Z;

        vertex = FINAL_tr.get_vertex_coord(i);
        FINAL3D(i+1,1) = vertex.X;  FINAL3D(i+1,2) = vertex.Y;  FINAL3D(i+1,3) = vertex.Z;
    }

    ORIG2D = ORIG3D * TRANS;
    FINAL2D = FINAL3D * TRANS2;

    W = triangle_strain(ORIG2D, FINAL2D, mu, kappa, indexSTRAINS, k_exp);

    return W;
}

void multivariate_histogram_normalization(MISCMATHS::BFMatrix& IN, MISCMATHS::BFMatrix& REF,
                                          const std::shared_ptr<newresampler::Mesh>& EXCL_IN,
                                          const std::shared_ptr<newresampler::Mesh>& EXCL_REF,
                                          int nthreads) {

    const int execution_threads = nthreads > IN.Nrows() ? IN.Nrows() : nthreads;

    #pragma omp parallel for num_threads(execution_threads)
    for(int d = 1; d <= (int) IN.Nrows(); d++)
    {
        NEWMAT::ColumnVector datain(IN.Ncols()), dataref(REF.Ncols());
        double max, min;
        NEWMAT::ColumnVector excluded_in(IN.Ncols()); excluded_in = 1;
        NEWMAT::ColumnVector excluded_ref(REF.Ncols()); excluded_ref = 1;
        const int numbins = 256;

        for (unsigned int i = 1; i <= IN.Ncols(); i++)
        {
            datain(i) = IN.Peek(d,i);

            if(EXCL_IN)
            { // if using an exclusion mask these values will be eliminated from the histogram matching
                if(EXCL_IN->get_dimension() >= d)
                    excluded_in(i) = EXCL_IN->get_pvalue(i-1,d-1);
                else
                    excluded_in(i) = EXCL_IN->get_pvalue(i-1);
            }
        }

        for (unsigned int i = 1; i <= REF.Ncols(); i++)
        {
            dataref(i) = REF.Peek(d,i);

            if(EXCL_REF)
            {
                if(EXCL_REF->get_dimension() >= d)
                    excluded_ref(i) = EXCL_REF->get_pvalue(i-1,d-1);
                else
                    excluded_ref(i) = EXCL_REF->get_pvalue(i-1);
            }
        }

        MISCMATHS::Histogram Hist_in(datain,numbins),
                             Hist_ref(dataref,numbins);
        Hist_in.generate(excluded_in);
        Hist_in.setexclusion(excluded_in);

        Hist_ref.generate(excluded_ref);
        Hist_ref.setexclusion(excluded_ref);

        Hist_in.generateCDF();
        Hist_ref.generateCDF();

        Hist_in.match(Hist_ref);

        datain = Hist_in.getsourceData();

        for (unsigned int i = 1; i <= IN.Ncols(); i++)
            IN.Set(d, i, datain(i));
    }
}

void set_data(const std::string& dataname, std::shared_ptr<MISCMATHS::BFMatrix>& BF, newresampler::Mesh& M, bool issparse) {

    if(issparse)
    {
        MISCMATHS::SpMat<double> sparse_mat(dataname);
        BF = std::shared_ptr<MISCMATHS::BFMatrix>(new MISCMATHS::SparseBFMatrix<double>(sparse_mat));
        // may not be desirable to do it this way if dimensions are v high
    }
    else
    {
        M.load(dataname,false,false);
        BF = std::shared_ptr<MISCMATHS::BFMatrix>(new MISCMATHS::FullBFMatrix(M.get_pvalues()));
    }

    if((int)BF->Ncols() != M.npvalues())
    {
        if((int)BF->Nrows() != M.npvalues())
            throw newresampler::MeshException("data does not match mesh dimensions");
        else
            BF = BF->Transpose();
    }
}

template<typename Iterator> inline bool next_combination(const Iterator first, Iterator k, const Iterator last) {
    if ((first == last) || (first == k) || (last == k))
        return false;
    Iterator itr1 = first;
    Iterator itr2 = last;
    ++itr1;
    if (last == itr1)
        return false;
    itr1 = last;
    --itr1;
    itr1 = k;
    --itr2;
    while (first != itr1)
    {
        if (*--itr1 < *itr2)
        {
            Iterator j = k;
            while (*itr1 >= *j) ++j;
            std::iter_swap(itr1,j);
            ++itr1;
            ++j;
            itr2 = k;
            std::rotate(itr1,j,last);
            while (last != j)
            {
                ++j;
                ++itr2;
            }
            std::rotate(k,itr2,last);
            return true;
        }
    }
    std::rotate(first,k,last);
    return false;
}

} //namespace newmeshreg
