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
#ifndef NEWRESAMPLER_MPOINT_H
#define NEWRESAMPLER_MPOINT_H

#include <vector>
#include <algorithm>

#include "point.h"

namespace newresampler {

class Mpoint {

    Point coord{};         // vertex coordinate
    int no{};              // counting index of vertex (where it appears in the points array of mesh)
    std::vector<int> trID; // faces adjacent to vertex
    std::vector<int> nID;  // neighbours of this vertex

public:
    Mpoint() = default;
    Mpoint(double x, double y, double z, int counter) : coord(Point(x, y, z)), no(counter) {}
    Mpoint(const Point p, int counter) : coord(p), no(counter) {}

    //---ACCESS---//
    const Point& get_coord() const { return coord; }
    int get_no() const { return no; }
    int ntriangles() const { return (int) trID.size(); }
    int nneighbours() const { return (int) nID.size(); }
    int get_trID(int i) const { return trID[i]; }
    const std::vector<int>& get_trID() const { return trID; }
    int get_nID(int i) const { return nID[i]; }
    const std::vector<int>& get_nID() const { return nID; }
    bool is_triangle(int j);
    bool is_neighbour(int j);

    //---MANIPULATION---//
    void set_coord(const Point& pt) { coord = pt; }
    void push_triangle(int i) { trID.push_back(i); }
    void push_neighbour(int i) { nID.push_back(i); }
    void normalize() { coord.normalize(); }
    void clear(){ trID.clear(); nID.clear(); }
};

//---Mpoint operators---//
bool operator==(const Mpoint &p1, const Mpoint &p2);
bool operator!=(const Mpoint &p1, const Mpoint &p2);
bool operator==(const Mpoint &p2, const Point &p1);
Point operator-(const Mpoint &p1, const Mpoint &p2);
Point operator-(const Point &p1, const Mpoint &p2);

} //namespace newresampler

#endif //NEWRESAMPLER_MPOINT_H
