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
#include "triangle.h"

namespace newresampler {

Triangle::Triangle(Point p1, Point p2, Point p3, int no) : no(no) {
    //TODO why all 0?
    vertices.emplace_back(std::make_shared<Mpoint>(p1, 0));
    vertices.emplace_back(std::make_shared<Mpoint>(p2, 0));
    vertices.emplace_back(std::make_shared<Mpoint>(p3, 0));
    area = calc_area();
}

Triangle::Triangle(const std::shared_ptr<Mpoint> &p1, const std::shared_ptr<Mpoint> &p2,
                   const std::shared_ptr<Mpoint> &p3, int no) : no(no) {
    vertices.push_back(p1);
    vertices.push_back(p2);
    vertices.push_back(p3);
    area = calc_area();
}

Point Triangle::normal() const {
    Point result = (vertices[2]->get_coord() - vertices[0]->get_coord()) *
                   (vertices[1]->get_coord() - vertices[0]->get_coord());
    result.normalize();
    return result;
}

double Triangle::calc_area() const {
    return 0.5 * ((vertices[2]->get_coord() - vertices[0]->get_coord()) *
                 (vertices[1]->get_coord() - vertices[0]->get_coord())).norm();
}

void Triangle::set(Point p1, Point p2, Point p3, int number)
{
    vertices.clear();
    vertices.emplace_back(std::make_shared<Mpoint>(p1,0));
    vertices.emplace_back(std::make_shared<Mpoint>(p2,1));
    vertices.emplace_back(std::make_shared<Mpoint>(p3,2));
    no = number;
    area = calc_area();
}

void Triangle::set_vertex(int i, const Point& p) {
    //TODO why 0?
    if(vertices.empty()) vertices.resize(3);
    vertices[i] = std::make_shared<Mpoint>(p,0);
}

std::vector<double> Triangle::get_angles() const { // get angles in order of vertex: 0,1,2

    std::vector<double> face_angles;

    Point v0 = vertices[2]->get_coord() - vertices[0]->get_coord();  // edge from vertex 0 to 2
    Point v1 = vertices[1]->get_coord() - vertices[0]->get_coord();  // edge from vertex 0 to 1
    Point v2 = vertices[2]->get_coord() - vertices[1]->get_coord();  // edge from vertex 1 to 2

    face_angles.emplace_back(acos((v0 | v1) / (v0.norm() * v1.norm())));
    face_angles.emplace_back(acos((v0 | v2) / (v0.norm() * v2.norm())));
    face_angles.emplace_back(acos((v1 | v2) / (v2.norm() * v1.norm())));

    return face_angles;
}

double Triangle::dist_to_point(const Point& x0) const {

    const Point x1(vertices[0]->get_coord());
    const Point x2(vertices[1]->get_coord());
    const Point x3(vertices[2]->get_coord());

    Point u;
    double d = 0.0;

    double dmin = std::numeric_limits<double>::max();
    // test edges
    u = x2 - x1;
    if (((x0 - x1) | u) > 0 && ((x0 - x2) | u) < 0) {
        d = (((x0 - x1) * (x0 - x2)).norm() / (x2 - x1).norm());
        if (d < dmin) dmin = d;
    }

    u = x3 - x1;
    if (((x0 - x1) | u) > 0 && ((x0 - x3) | u) < 0) {
        d = (((x0 - x1) * (x0 - x3)).norm() / (x3 - x1).norm());
        if (d < dmin) dmin = d;
    }

    u = x3 - x2;
    if (((x0 - x2) | u) > 0 && ((x0 - x3) | u) < 0) {
        d = (((x0 - x2) * (x0 - x3)).norm() / (x3 - x2).norm());
        if (d < dmin) dmin = d;
    }

    d = (x0 - x1).norm();
    if (d < dmin) dmin = d;
    d = (x0 - x2).norm();
    if (d < dmin) dmin = d;
    d = (x0 - x3).norm();
    if (d < dmin) dmin = d;

    return dmin;
}

std::map<int,double> calc_barycentric_weights(const Point& v1, const Point& v2,
                                              const Point& v3, const Point& vref,
                                              int n1, int n2, int n3){

    std::map<int,double> weights;

    Point PP = project_point(vref, v1, v2, v3);

    double Aa = compute_area(PP, v2, v3);
    double Ab = compute_area(PP, v1, v3);
    double Ac = compute_area(PP, v1, v2);

    double A = Aa + Ab + Ac;

    weights[n1] = Aa / A;
    weights[n2] = Ab / A;
    weights[n3] = Ac / A;

    return weights;
}

double barycentric_interpolation(const Point& v1, const Point& v2, const Point& v3, const Point& vref, double va1, double va2, double va3) {

    double Aa = compute_area(vref, v2, v3);
    double Ab = compute_area(vref, v1, v3);
    double Ac = compute_area(vref, v1, v2);

    double A = Aa + Ab + Ac;
    Aa = Aa / A;
    Ab = Ab / A;
    Ac = Ac / A;

    return Aa * va1 + Ab * va2 + Ac * va3;
}

Point barycentric(const Point& v1, const Point& v2, const Point& v3, const Point& vref, const Point& va1, const Point& va2, const Point& va3){
    //TODO rename and check

    double Aa = compute_area(vref, v2, v3);
    double Ab = compute_area(vref, v1, v3);
    double Ac = compute_area(vref, v1, v2);

    double A = Aa + Ab + Ac;
    Aa = Aa / A;
    Ab = Ab / A;
    Ac = Ac / A;

    return  va1 * Aa + va2 * Ab + va3 * Ac;
}

} //namespace newresampler
