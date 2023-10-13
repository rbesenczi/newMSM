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

Octree search code used with the permission of Tim Coalson under the same licence as above.
Copyright (C) Washington University School of Medicine
Original: https://github.com/Washington-University/workbench/blob/master/src/Files/SignedDistanceHelper.h
https://github.com/Washington-University/workbench/blob/master/src/Files/SignedDistanceHelper.cxx
*/
#ifndef NEWRESAMPLER_OCTREE_H
#define NEWRESAMPLER_OCTREE_H

#include "node.h"
#include "mesh.h"

namespace newresampler {

class Octree {

    Node* octree_root;

    void initialize_tree(const std::vector<Triangle>& triangles);
    void add_triangle(Node* node, const Triangle& triangle, const Point& lower_bounds, const Point& upper_bounds);
    double distance_to_triangle(const Point& pt, const Triangle& tr) const;

public:
    explicit Octree(const Mesh& target);
    ~Octree() { delete octree_root; octree_root = nullptr; }

    Triangle get_closest_triangle(const Point& pt) const;
    int get_closest_vertex_ID(const Point& pt) const;

    //disable copy and move for now
    Octree(const Octree&) = delete;
    Octree& operator=(const Octree&) = delete;
    Octree(Octree&&) = delete;
    Octree& operator=(Octree&&) = delete;
};

} // namespace newresampler

#endif //NEWRESAMPLER_OCTREE_H
