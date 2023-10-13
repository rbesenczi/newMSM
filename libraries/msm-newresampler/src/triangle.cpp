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
    vertices.emplace_back(std::make_shared<Mpoint>(p1, 0));
    vertices.emplace_back(std::make_shared<Mpoint>(p2, 0));
    vertices.emplace_back(std::make_shared<Mpoint>(p3, 0));
    area = calc_area();
}

Triangle::Triangle(const std::shared_ptr<Mpoint> &p1, const std::shared_ptr<Mpoint> &p2,
                   const std::shared_ptr<Mpoint> &p3, int number) {
    vertices.push_back(p1);
    vertices.push_back(p2);
    vertices.push_back(p3);
    no = number;
    area = calc_area();
}

Triangle Triangle::copy() const {
    Triangle t;
    for (const auto& i: vertices)
        t.vertices.emplace_back(std::make_shared<Mpoint>(*i));
    t.no = no;
    t.area = area;
    return t;
}

Point Triangle::centroid() const {
    return {(vertices[0]->get_coord().X + vertices[1]->get_coord().X + vertices[2]->get_coord().X) / 3,
            (vertices[0]->get_coord().Y + vertices[1]->get_coord().Y + vertices[2]->get_coord().Y) / 3,
            (vertices[0]->get_coord().Z + vertices[1]->get_coord().Z + vertices[2]->get_coord().Z) / 3};
}

Point Triangle::normal() const {
    Point result = (vertices[2]->get_coord() - vertices[0]->get_coord()) *
                   (vertices[1]->get_coord() - vertices[0]->get_coord());
    result.normalize();
    return result;
}

double Triangle::calc_area() const {
    Point result = (vertices[2]->get_coord() - vertices[0]->get_coord()) *
                   (vertices[1]->get_coord() - vertices[0]->get_coord());
    return 0.5 * result.norm();
}

void Triangle::set(Point p1, Point p2, Point p3, int index)
{
    std::shared_ptr<Mpoint> m1 = std::make_shared<Mpoint>(p1,0);
    std::shared_ptr<Mpoint> m2 = std::make_shared<Mpoint>(p2,1);
    std::shared_ptr<Mpoint> m3 = std::make_shared<Mpoint>(p3,2);
    vertices.clear(); vertices.push_back(m1); vertices.push_back(m2); vertices.push_back(m3);
    no = index;
    area = calc_area();
}

void Triangle::set_vertex(int i, const Point& p) {
    std::shared_ptr<Mpoint> m = std::make_shared<Mpoint>(p,0);
    if(vertices.empty()) vertices.resize(3);
    vertices[i] = m;
}

std::vector<double> Triangle::get_angles() const { // get angles in order of vertex: 0,1,2

    std::vector<double> face_angles;

    Point v0 = vertices[2]->get_coord() - vertices[0]->get_coord();  // edge from vertex 0 to 2
    Point v1 = vertices[1]->get_coord() - vertices[0]->get_coord();  // edge from vertex 0 to 1
    Point v2 = vertices[2]->get_coord() - vertices[1]->get_coord();  // edge from vertex 1 to 2

    double dot01 = v0 | v1;
    double dot02 = v0 | v2;
    double dot12 = v1 | v2;

    face_angles.emplace_back(acos(dot01 / (v0.norm() * v1.norm())));
    face_angles.emplace_back(acos(dot02 / (v0.norm() * v2.norm())));
    face_angles.emplace_back(acos(dot12 / (v2.norm() * v1.norm())));

    return face_angles;
}

bool Triangle::is_inside(const Point& x) const {

    Point v0 = vertices[2]->get_coord() - vertices[0]->get_coord();
    Point v1 = vertices[1]->get_coord() - vertices[0]->get_coord();
    Point v2 = x - vertices[0]->get_coord();

    double dot00 = v0 | v0;
    double dot01 = v0 | v1;
    double dot02 = v0 | v2;
    double dot11 = v1 | v1;
    double dot12 = v1 | v2;

    double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);

    double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    return (u > 0) && (v > 0) && (u + v < 1);
}

double Triangle::dist_to_point(const Point& x0) const {

    double d = 0.0;

    Point x1(vertices[0]->get_coord());
    Point x2(vertices[1]->get_coord());
    Point x3(vertices[2]->get_coord());

    Point u;

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

    double A, Aa, Ab, Ac;
    Point PP;
    std::map<int,double> weights;

    project_point(vref, v1, v2, v3, PP);

    Aa = compute_area(PP, v2, v3);
    Ab = compute_area(PP, v1, v3);
    Ac = compute_area(PP, v1, v2);

    A = Aa + Ab + Ac;

    weights[n1] = Aa / A;
    weights[n2] = Ab / A;
    weights[n3] = Ac / A;

    return weights;
}

double barycentric_weight(const Point& v1, const Point& v2, const Point& v3, const Point& vref, double va1, double va2, double va3) {

    double Aa = compute_area(vref, v2, v3);
    double Ab = compute_area(vref, v1, v3);
    double Ac = compute_area(vref, v1, v2);

    double A = Aa + Ab + Ac;
    Aa = Aa / A;
    Ab = Ab / A;
    Ac = Ac / A;

    return  Aa * va1 + Ab * va2 + Ac * va3;
}

Point barycentric(const Point& v1, const Point& v2, const Point& v3, const Point& vref, const Point& va1, const Point& va2, const Point& va3){

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
