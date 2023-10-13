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
*/
#include "point.h"

namespace newresampler {

void Point::normalize() {
    double n = norm();
    if (n != 0)
    {
        X /= n;
        Y /= n;
        Z /= n;
    }
}

bool same_side(const Point &p1, const Point &p2, const Point &a, const Point &b) {

    Point cp1 = (b - a) * (p1 - a);
    Point cp2 = (b - a) * (p2 - a);

    if ((cp1 | cp2) >= -1E-6) return true; else return false;
}

bool point_in_triangle(const Point &p, const Point &a, const Point &b, const Point &c) {

    if (same_side(p, a, b, c) && same_side(p, b, c, a) && same_side(p, c, a, b))
        return true;
    else return false;
}

Point project_point(const Point &vb, const Point &v1, const Point &v2,
                    const Point &v3, Point &PP) {

    Point s1 = v3 - v1;
    s1.normalize();
    Point s2 = v2 - v1;
    s2.normalize();
    Point s3 = s1 * s2; // s3= normal
    s3.normalize();

    double si = (s3 | v1) / (s3 | vb);
    // formula for line plane intersection s3.(v1-sivb)=0
    // i.e. (v1-sivb) should be perpendicular to plane normal s3

    PP = vb * si;  // PP therefore where vb intersects plane
    return s3;
}

void project_point(const Point& v1, const Point& v2, const Point& v3, Point& PP)
{
    Point s1 = v3 - v1; s1.normalize();
    Point s2 = v2 - v1; s2.normalize();
    Point s3 = s1 * s2; s3.normalize();

    double si = (s3|v1)/(s3|PP); // formula for line plane intersection s3.(v1-sivb)=0 i.e. (v1-sivb) should be perpendicular to plane normal s3

    PP *= si;  // PP therefore where vb intersects plane
}

void project_point(const Point& vb, const Tangs& T, double& e1coord, double& e2coord) {
    e1coord = vb | T.e1;
    e2coord = vb | T.e2;
}

double compute_area(const Point &v0, const Point &v1, const Point &v2) {

    Point res, v0v1, v0v2;

    v0v1 = v1 - v0;
    v0v2 = v2 - v0;
    res = v0v1 * v0v2;

    return 0.5 * res.norm();
}

NEWMAT::Matrix form_matrix_from_points(const Point& p1, const Point& p2, const Point& p3, bool trans) {

    NEWMAT::Matrix T(3,3);

    T(1,1) = p1.X;
    T(1,2) = p2.X;
    T(1,3) = p3.X;
    T(2,1) = p1.Y;
    T(2,2) = p2.Y;
    T(2,3) = p3.Y;
    T(3,1) = p1.Z;
    T(3,2) = p2.Z;
    T(3,3) = p3.Z;

    if(trans)
        return T.t();
    else
        return T;
}

NEWMAT::Matrix estimate_rotation_matrix(const Point& p1, const Point& p2) {
    // ci is start point index is end point
    NEWMAT::Matrix I(3,3),u(3,3),R(3,3);
    Point ci = p1;
    Point index = p2;
    ci.normalize();
    index.normalize();
    double theta = acos(ci | index);
    const double PI = 3.14159265359, EPSILON = 1.0E-8;

    if (theta > PI)
        throw MeshException("rotation angle is greater than 90 degrees");

    I << 1 << 0 << 0
      << 0 << 1 << 0
      << 0 << 0 << 1;

    //calculate axis of rotation from cross product between centre point and neighbour - could also use barycentres
    Point cross = ci * index;
    cross.normalize();

    if (fabs(1 - (ci | index)) < EPSILON)
        R = I;
    else if (cross.norm() == 0)
    {
        R << -1 << 0 << 0
          << 0 << -1 << 0
          << 0 << 0 << -1;
    }
    else
    {
        u << 0 << -cross.Z << cross.Y
          << cross.Z << 0 << -cross.X
          << -cross.Y << cross.X << 0;

        // now use rodriguez formula to rotate all data points that contribute the similarity of this control point
        if (fabs(-1 - (ci | index)) < EPSILON) { // numerical problems at theta cose to PI
            NEWMAT::Matrix outerprod(3, 3);

            outerprod << cross.X * cross.X << cross.X * cross.Y << cross.X * cross.Z
                      << cross.Y * cross.X << cross.Y * cross.Y << cross.Y * cross.Z
                      << cross.Z * cross.X << cross.Z * cross.Y << cross.Z * cross.Z;

            R = 2 * outerprod - I;
        }
        else
            R = I + u * sin(theta) + (1 - cos(theta)) * (u * u);
    }

    return R;
}

NEWMAT::ReturnMatrix euler_rotate(const NEWMAT::ColumnVector& vector, double w1, double w2, double w3) {

    NEWMAT::Matrix rotation(3, 3);
    rotation
            << cos(w2) * cos(w3)
            << -cos(w1) * sin(w3) + sin(w1) * sin(w2) * cos(w3)
            << sin(w1) * sin(w3) + cos(w1) * sin(w2) * cos(w3)
            << cos(w2) * sin(w3)
            << cos(w1) * cos(w3) + sin(w1) * sin(w2) * sin(w3)
            << -sin(w1) * cos(w3) + cos(w1) * sin(w2) * sin(w3)
            << -sin(w2)
            << sin(w1) * cos(w2) << cos(w1) * cos(w2);

    NEWMAT::ColumnVector vector_rot = rotation.t() * vector;

    vector_rot.Release();
    return vector_rot;
}

double operator|(const Point &v1, const Point &v2) {
    return v1.X * v2.X + v1.Y * v2.Y + v1.Z * v2.Z;
}

Point operator*(const Point &v1, const Point &v2) {

    return {v1.Y * v2.Z - v1.Z * v2.Y,
            v2.X * v1.Z - v2.Z * v1.X,
            v1.X * v2.Y - v2.X * v1.Y};
}

Point operator+(const Point &v1, const Point &v2) {
    return {v1.X + v2.X, v1.Y + v2.Y, v1.Z + v2.Z};
}

Point operator-(const Point &v1, const Point &v2) {
    return {v1.X - v2.X, v1.Y - v2.Y, v1.Z - v2.Z};
}

Point operator/(const Point &v, const double &d) {
    if (d != 0)
        return {v.X / d, v.Y / d, v.Z / d};
    else
        throw MeshException{"newresampler::Point operator/ division by zero"};
}

Point operator*(const Point &v, const double &d) {
    return {v.X * d, v.Y * d, v.Z * d};
}

Point operator*(const NEWMAT::Matrix &M, const Point &v) {
    if (M.Ncols() != 3 || M.Nrows() != 3)
        throw MeshException{"newresampler::Point matrix multiply error: matrix should be 3x3"};
    return {M(1, 1) * v.X + M(1, 2) * v.Y + M(1, 3) * v.Z,
            M(2, 1) * v.X + M(2, 2) * v.Y + M(2, 3) * v.Z,
            M(3, 1) * v.X + M(3, 2) * v.Y + M(3, 3) * v.Z};
}

Point operator*(const Point &v, const NEWMAT::Matrix &M) {
    if (M.Ncols() != 3 || M.Nrows() != 3)
        throw MeshException{"newresampler::Point matrix multiply error: matrix should be 3x3"};
    return {v.X * M(1, 1) + v.Y * M(2, 1) + v.Z * M(3, 1),
            v.X * M(1, 2) + v.Y * M(2, 2) + v.Z * M(3, 2),
            v.X * M(1, 3) + v.Y * M(2, 3) + v.Z * M(3, 3)};
}

Point operator+=(Point &p1, const Point &p2) {
    p1.X += p2.X;
    p1.Y += p2.Y;
    p1.Z += p2.Z;
    return p1;
}

Point operator*=(Point &p, const double &d) {
    p.X *= d;
    p.Y *= d;
    p.Z *= d;
    return p;
}

Point operator/=(Point &p, const double &d) {
    if (d != 0) {
        p.X /= d;
        p.Y /= d;
        p.Z /= d;
    } else throw MeshException{"newresampler::Point operator/= division by zero"};
    return p;
}

bool operator==(const Point &p1, const Point &p2) {
    return ((std::fabs(p1.X - p2.X) < 1e-4)
        && (std::fabs(p1.Y - p2.Y) < 1e-4)
        && (std::fabs(p1.Z - p2.Z) < 1e-4));
}

bool operator!=(const Point& p1, const Point& p2) {
    return !(p1==p2);
}

std::ostream& operator<<(std::ostream& os, const Point& pt) {
    return os << pt.X << ' ' << pt.Y << ' ' << pt.Z;
}

} //namespace newresampler
