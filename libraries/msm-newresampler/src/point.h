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
#ifndef NEWRESAMPLER_POINT_H
#define NEWRESAMPLER_POINT_H

#include <cmath>

#include "armawrap/newmat.h"

#include "meshException.h"

namespace newresampler {

class Point {

public:
    double X, Y, Z;

    Point(): X(0.0), Y(0.0), Z(0.0) {};
    Point(double x, double y, double z): X(x), Y(y), Z(z) {};

    inline double norm() const { return std::sqrt(X * X + Y * Y + Z * Z); }
    void normalize();
};

struct Tangs { Point e1, e2; };

//---Utility functions related with Point---//
bool same_side(const Point &, const Point &, const Point &, const Point &);
bool point_in_triangle(const Point &, const Point &, const Point &, const Point &);
Point project_point(const Point &, const Point &, const Point &, const Point &, Point &);
void project_point(const Point &, const Point &, const Point &, Point&);
void project_point(const Point& vb, const Tangs& T, double& e1coord, double& e2coord);
double compute_area(const Point &, const Point &, const Point &);
NEWMAT::Matrix form_matrix_from_points(const Point& p1, const Point& p2, const Point& p3, bool trans = false);
NEWMAT::Matrix estimate_rotation_matrix(const Point& p1, const Point& p2);
NEWMAT::ReturnMatrix euler_rotate(const NEWMAT::ColumnVector& vector, double w1, double w2, double w3);

//---Point operators---//
double operator|(const Point &p1, const Point &p2); // dot product
Point operator*(const Point &p1, const Point &p2); // cross product
Point operator*(const Point &v, const double &d);
Point operator/(const Point &v, const double &d);
Point operator-(const Point &p1, const Point &p2);
Point operator+(const Point &p1, const Point &p2);
Point operator*(const NEWMAT::Matrix &M, const Point &p1);
Point operator*(const Point &p1, const NEWMAT::Matrix &M);
Point operator+=(Point &p1, const Point &p2);
Point operator*=(Point &p, const double &d);
Point operator/=(Point &p, const double &d);
bool operator==(const Point &p1, const Point &p2);
bool operator!=(const Point &p1, const Point &p2);

std::ostream& operator<<(std::ostream&, const Point&);

} //namespace newresampler

#endif //NEWRESAMPLER_POINT_H
