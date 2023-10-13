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
#ifndef NEWRESAMPLER_TRIANGLE_H
#define NEWRESAMPLER_TRIANGLE_H

#include <memory>

#include "mpoint.h"

namespace newresampler {

class Triangle {

    std::vector<std::shared_ptr<Mpoint>> vertices;
    int no = -1;
    double area = 0.0;

public:
    Triangle(){};

    Triangle(Point p1, Point p2, Point p3, int no);
    Triangle(const std::shared_ptr<Mpoint>& p1,
             const std::shared_ptr<Mpoint>& p2,
             const std::shared_ptr<Mpoint>& p3,
             int number);

    //---ACCESS---//
    int get_no() const { return no; }
    double get_area() const { return area; }
    const Mpoint& get_vertex(int i) const { return *vertices[i]; }
    const std::shared_ptr<Mpoint>& get_vertex_ptr(int i) const { return vertices[i]; }
    const Point& get_vertex_coord(int i) const { return vertices[i]->get_coord(); }
    int get_vertex_no(int i) const { return vertices[i]->get_no(); }

    //---MANIPULATION---//
    inline void swap_orientation() { std::swap(vertices[1], vertices[2]); }
    void set(Point p1, Point p2, Point p3, int index);
    void set_vertex(int i, const Point& p);

    //---UTILITY---//
    Triangle copy() const;
    Point centroid() const;
    Point normal() const;
    double calc_area() const;
    std::vector<double> get_angles() const;
    bool is_inside(const Point &x) const;
    double dist_to_point(const Point &x0) const;
};

std::map<int,double> calc_barycentric_weights(const Point&, const Point&, const Point&, const Point&, int, int, int);
double barycentric_weight(const Point& v1, const Point& v2, const Point& v3, const Point& vref, double va1, double va2, double va3);
Point barycentric(const Point& v1, const Point& v2, const Point& v3, const Point& vref, const Point& va1, const Point& va2, const Point& va3);

} //namespace newresampler

#endif //NEWRESAMPLER_TRIANGLE_H
