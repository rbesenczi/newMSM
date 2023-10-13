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
#include "octree.h"

namespace newresampler {

Octree::Octree(const Mesh& target) {

    const int limit = 101;
    double lower_bounds[3] = {-limit, -limit, -limit };
    double upper_bounds[3] = {limit, limit, limit };

    octree_root = new Node(lower_bounds, upper_bounds);

    initialize_tree(target.get_all_triangles());
}

void Octree::initialize_tree(const std::vector<Triangle>& triangles) {

    for(const auto& tr: triangles)
    {
        Point temp_lower_bound = tr.get_vertex_coord(0);
        Point temp_upper_bound = tr.get_vertex_coord(0);

        for (int vertex_id = 1; vertex_id < 3; ++vertex_id)
        {
            Point triangle_vertex = tr.get_vertex_coord(vertex_id);
            if (triangle_vertex.X < temp_lower_bound.X) temp_lower_bound.X = triangle_vertex.X;
            if (triangle_vertex.Y < temp_lower_bound.Y) temp_lower_bound.Y = triangle_vertex.Y;
            if (triangle_vertex.Z < temp_lower_bound.Z) temp_lower_bound.Z = triangle_vertex.Z;
            if (triangle_vertex.X > temp_upper_bound.X) temp_upper_bound.X = triangle_vertex.X;
            if (triangle_vertex.Y > temp_upper_bound.Y) temp_upper_bound.Y = triangle_vertex.Y;
            if (triangle_vertex.Z > temp_upper_bound.Z) temp_upper_bound.Z = triangle_vertex.Z;
        }
        add_triangle(octree_root, tr, temp_lower_bound, temp_upper_bound);
    }
}

void Octree::add_triangle(Node* node, const Triangle& triangle, const Point& lower_bounds, const Point& upper_bounds) {

    if (node->is_leaf)
    {
        node->push_triangle(triangle);
        int num_triangles = node->triangles_size();
        if (num_triangles >= MAX_TRIANGLES)
        {
            Point temp_lower_bounds, temp_upper_bounds;
            int total_size = 0;
            int num_split = 0;
            for (int i = 0; i < num_triangles; ++i)
            {
                Triangle temp_triangle = node->get_triangle(i);
                temp_lower_bounds = temp_triangle.get_vertex_coord(0);
                temp_upper_bounds = temp_triangle.get_vertex_coord(0);

                for (int vertex_id = 1; vertex_id < 3; ++vertex_id)
                {
                    Point temp_triangle_vertex = temp_triangle.get_vertex_coord(vertex_id);
                    if (temp_triangle_vertex.X < temp_lower_bounds.X) temp_lower_bounds.X = temp_triangle_vertex.X;
                    if (temp_triangle_vertex.Y < temp_lower_bounds.Y) temp_lower_bounds.Y = temp_triangle_vertex.Y;
                    if (temp_triangle_vertex.Z < temp_lower_bounds.Z) temp_lower_bounds.Z = temp_triangle_vertex.Z;
                    if (temp_triangle_vertex.X > temp_upper_bounds.X) temp_upper_bounds.X = temp_triangle_vertex.X;
                    if (temp_triangle_vertex.Y > temp_upper_bounds.Y) temp_upper_bounds.Y = temp_triangle_vertex.Y;
                    if (temp_triangle_vertex.Z > temp_upper_bounds.Z) temp_upper_bounds.Z = temp_triangle_vertex.Z;
                }

                std::vector<bool> min_oct = node->containing_oct(temp_lower_bounds);
                std::vector<bool> max_oct = node->containing_oct(temp_upper_bounds);

                int split_size = 8;

                for (int d = 0; d < 3; d++)
                    if (min_oct[d] == max_oct[d]) split_size >>= 1;

                total_size += split_size;
                if (split_size != 8) ++num_split;
            }
            if (num_split > 0 && total_size < 3 * (int)num_triangles)
            {
                node->make_children();
                for (int d = 0; d < num_triangles; ++d)
                {
                    Triangle temp_triangle = node->get_triangle(d);

                    temp_lower_bounds = temp_triangle.get_vertex_coord(0);
                    temp_upper_bounds = temp_triangle.get_vertex_coord(0);

                    for (int vertex_id = 1; vertex_id < 3; ++vertex_id)
                    {
                        Point temp_triangle_vertex = temp_triangle.get_vertex_coord(vertex_id);
                        if (temp_triangle_vertex.X < temp_lower_bounds.X) temp_lower_bounds.X = temp_triangle_vertex.X;
                        if (temp_triangle_vertex.Y < temp_lower_bounds.Y) temp_lower_bounds.Y = temp_triangle_vertex.Y;
                        if (temp_triangle_vertex.Z < temp_lower_bounds.Z) temp_lower_bounds.Z = temp_triangle_vertex.Z;
                        if (temp_triangle_vertex.X > temp_upper_bounds.X) temp_upper_bounds.X = temp_triangle_vertex.X;
                        if (temp_triangle_vertex.Y > temp_upper_bounds.Y) temp_upper_bounds.Y = temp_triangle_vertex.Y;
                        if (temp_triangle_vertex.Z > temp_upper_bounds.Z) temp_upper_bounds.Z = temp_triangle_vertex.Z;
                    }

                    for (auto& i: node->children)
                        for (auto& j: i)
                            for (auto& k: j)
                                if (k->can_contain(temp_lower_bounds, temp_upper_bounds))
                                    add_triangle(k, temp_triangle, temp_lower_bounds, temp_upper_bounds);
                }
                node->clear_triangles();
            }
        }
    }
    else
    {
        for (auto& i: node->children)
            for (auto& j: i)
                for (auto& k: j)
                    if (k->can_contain(lower_bounds, upper_bounds))
                        add_triangle(k, triangle, lower_bounds, upper_bounds);
    }
}

double Octree::distance_to_triangle(const Point& pt, const Triangle& tr) const {

    Point mP;
    const Point& v0 = tr.get_vertex_coord(0),
                 v1 = tr.get_vertex_coord(1),
                 v2 = tr.get_vertex_coord(2);

    project_point(pt, v0, v1, v2, mP);

    if (point_in_triangle(mP, v0, v1, v2))
        return tr.dist_to_point(mP);
    else return -1.0;   //don't like magic numbers, but this indicates that the point is not in the triangle
}

Triangle Octree::get_closest_triangle(const Point &pt) const {

    if(!octree_root->contains_point(pt)) { throw MeshException("Point is not in the bounding box of the mesh"); }

    Triangle closest_triangle;
    double best_distance = std::numeric_limits<double>::max();

    Node* current_oct = octree_root;

    while (!current_oct->is_leaf)
        for (auto& i: current_oct->children)
            for (auto& j: i)
                for (auto& k: j)
                    if (k->contains_point(pt))
                        current_oct = k;

    for (int i = 0; i < current_oct->triangles_size(); ++i) {
        double temp_distance = distance_to_triangle(pt, current_oct->get_triangle(i));
        if(temp_distance > -1.0 && temp_distance < best_distance) {
            closest_triangle = current_oct->get_triangle(i);
            best_distance = temp_distance;
        }
    }

    if(closest_triangle.get_no() == -1) {
        best_distance = std::numeric_limits<double>::max();
        for (auto &i: current_oct->parent->children)
            for (auto &j: i)
                for (auto &close_oct: j)
                    for (int tr = 0; tr < close_oct->triangles_size(); ++tr) {
                        double temp_distance = distance_to_triangle(pt, close_oct->get_triangle(tr));
                        if (temp_distance > -1.0 && temp_distance < best_distance) {
                            closest_triangle = close_oct->get_triangle(tr);
                            best_distance = temp_distance;
                        }
                    }
    }

    if(closest_triangle.get_no() == -1) {
        //This is an unfortunate case... Probably needs some more refinement.
        best_distance = std::numeric_limits<double>::max();
        constexpr double RAD = 100.0;
        for (auto &i: current_oct->parent->children)
            for (auto &j: i)
                for (auto close_oct: j)
                    for (int tr = 0; tr < close_oct->triangles_size(); ++tr)
                        for (int v = 0; v < 3; ++v) {
                            Point CP = close_oct->get_triangle(tr).get_vertex_coord(v);
                            double temp_distance = 2 * RAD * asin((CP - pt).norm() / (2 * RAD));
                            if (temp_distance < best_distance) {
                                closest_triangle = close_oct->get_triangle(tr);
                                best_distance = temp_distance;
                            }
                        }
    }

    if(closest_triangle.get_no() == -1) throw MeshException("Error in octree...");

    return closest_triangle;
}

int Octree::get_closest_vertex_ID(const Point &pt) const {

    Triangle closest_triangle = get_closest_triangle(pt);

    double dist = std::numeric_limits<double>::max();
    int closest_vertex_ID = 0;

    for(int v = 0; v < 3; ++v)
    {
        double current_dist = (pt - closest_triangle.get_vertex_coord(v)).norm();
        if(current_dist < dist)
        {
            closest_vertex_ID = closest_triangle.get_vertex_no(v);
            dist = current_dist;
        }
    }

    return closest_vertex_ID;
}

} //namespace newresampler
