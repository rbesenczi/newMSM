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
Original: https://github.com/Washington-University/workbench/blob/master/src/Common/OctTree.h
*/
#ifndef NEWRESAMPLER_NODE_H
#define NEWRESAMPLER_NODE_H

#include "triangle.h"

namespace newresampler {

static const int MAX_TRIANGLES = 50; //maximum number of triangles in one node of the tree

class Node {

    std::vector<Triangle> triangles; //data if this is a leaf, empty otherwise

public:
    Node* children[2][2][2]; //pointers to children
    Node* parent;
    double bounds[3][3]; //the bounds of this node
    bool is_leaf = true;

    Node();
    Node(const double lower_bounds[3], const double upper_bounds[3]);
    ~Node();

    int triangles_size() const { return (int) triangles.size(); }
    const auto& get_triangle(int i) const { return triangles[i]; }
    void push_triangle(const Triangle& tr) { triangles.push_back(tr); }
    void clear_triangles() { triangles.clear(); }

    void make_children();
    bool can_contain(const Point& lower_bounds, const Point& upper_bounds);
    std::vector<bool> containing_oct(const Point &pt);
    bool contains_point(const Point &pt);

    //disable copy and move for now
    Node(const Node&) = delete;
    Node& operator=(const Node&) = delete;
    Node(Node&&) = delete;
    Node& operator=(Node&&) = delete;
};

} //namespace newresampler

#endif //NEWRESAMPLER_NODE_H
