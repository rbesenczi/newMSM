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
#include "node.h"

namespace newresampler {

Node::Node() {

    for(auto& i : children)
        for(auto& j : i)
            for (auto& k: j)
                k = nullptr;

    triangles.reserve(MAX_TRIANGLES);
    is_leaf = true;
}

Node::Node(const double lower_bounds[3], const double upper_bounds[3]) {

    for(auto& i : children)
        for(auto& j : i)
            for (auto& k: j)
                k = nullptr;

    triangles.reserve(MAX_TRIANGLES);
    is_leaf = true;

    for (int i = 0; i < 3; ++i)
    {
        bounds[i][0] = lower_bounds[i];
        bounds[i][2] = upper_bounds[i];
        bounds[i][1] = (bounds[i][0] + bounds[i][2]) / 2.0;
    }
}

Node::~Node() {

    is_leaf = true;

    for(auto& i : children)
        for(auto& j : i)
            for(auto& k : j)
            {
                delete k; k = nullptr;
            }
}

bool Node::contains_point(const Point& pt) {

    double point[3] = { pt.X, pt.Y, pt.Z };

    for (int i = 0; i < 3; ++i)
    {
        if (point[i] < bounds[i][0]) return false;
        if (point[i] > bounds[i][2]) return false;
    }
    return true;
}

std::vector<bool> Node::containing_oct(const Point &pt) {

    std::vector<bool> child_octs(3);

    double point[3] = {pt.X,pt.Y,pt.Z };

    for (int i = 0; i < 3; ++i)
        child_octs[i] = point[i] < bounds[i][1];

    return child_octs;
}

void Node::make_children() {

    is_leaf = false;

    int child_oct[3];
    for (child_oct[0] = 0; child_oct[0] < 2; ++child_oct[0])
        for (child_oct[1] = 0; child_oct[1] < 2; ++child_oct[1])
            for (child_oct[2] = 0; child_oct[2] < 2; ++child_oct[2])
            {
                auto* temp = new Node();
                for (int i = 0; i < 3; ++i)
                {
                    temp->bounds[i][0] = bounds[i][child_oct[i]];
                    temp->bounds[i][2] = bounds[i][child_oct[i] + 1];
                    temp->bounds[i][1] = (temp->bounds[i][0] + temp->bounds[i][2]) / 2.0;
                }
                temp->parent = this;
                children[child_oct[0]][child_oct[1]][child_oct[2]] = temp;
            }
}

bool Node::can_contain(const Point& lower_bounds, const Point& upper_bounds) {

    double temp_min[3] = {lower_bounds.X, lower_bounds.Y, lower_bounds.Z};
    double temp_max[3] = {upper_bounds.X, upper_bounds.Y, upper_bounds.Z};

    for (int i = 0; i < 3; ++i)
        if (temp_max[i] < bounds[i][0] || temp_min[i] > bounds[i][2]) return false;
    return true;
}

} //namespace newresampler
