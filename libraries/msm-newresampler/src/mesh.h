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
#ifndef NEWRESAMPLER_MESH_H
#define NEWRESAMPLER_MESH_H

#include <cstdio>
#include <iomanip>
#include <omp.h>

#include "miscmaths/miscmaths.h"
#include "miscmaths/bfmatrix.h"
#include "newmesh/giftiInterface.h"

#include "triangle.h"

namespace newresampler {

class Mesh {

    std::vector<std::shared_ptr<Mpoint>> points;
    std::vector<Triangle> triangles;
    std::vector<Point> normals;

    std::vector<std::vector<float>> pvalues;
    std::vector<std::vector<float>> tvalues;

    //---METADATA---//
    std::vector<NEWMESH::GIFTImeta> global_metaData;
    std::vector<NEWMESH::GIFTImeta> global_Attributes;
    std::map<int, NEWMESH::GIFTIlabel> global_GIFTIlabels;
    std::vector<NEWMESH::GIFTIcoordinateSystem> global_defaultcoord;

public:
    Mesh();

    //---COPY/MOVE---//
    Mesh(const Mesh &m);
    Mesh& operator=(const Mesh& m);
    Mesh(Mesh&& m) noexcept;
    Mesh& operator=(Mesh&& m) noexcept;

    ~Mesh() = default;

    //---ACCESS---//
    int nvertices() const { return (int)points.size(); }
    int ntriangles() const { return (int)triangles.size(); }
    float get_pvalue(int i, int dim = 0) const;
    int get_dimension() const { return (int)pvalues.size(); }
    int npvalues (int dim = 0) const { return (int) pvalues[dim].size(); }
    const Point& get_coord(int n) const;
    const Mpoint& get_point(int n) const { return *points[n]; }
    const std::vector<Triangle>& get_all_triangles() const { return triangles; }
    int get_total_triangles(int i) const;
    double get_triangle_area(int tID) const { return triangles[tID].get_area(); }
    NEWMAT::Matrix get_pvalues() const;

    const Triangle& get_triangle(int n) const {
        if (n >= (int) triangles.size() || triangles.empty())
            throw MeshException("get_triangle: index exceeds face dimensions");
        return triangles[n];
    }

    const Triangle& get_triangle_from_vertex(int n, int ID) const {
        if (n >= (int) triangles.size() || (int) triangles.empty())
            throw MeshException("get_triangle: index exceeds face dimensions");
        return triangles[points[n]->get_trID(ID)];
    }

    const Point& get_triangle_vertex(int n, int i) const {
        if (n >= (int) triangles.size() || triangles.empty())
            throw MeshException("get_triangle: index exceeds face dimensions");
        return triangles[n].get_vertex_coord(i);
    }

    const Point& get_normal(int n) const {
        if (normals.size() < points.size())
            throw MeshException("get_normal: normals have not been calculated, apply estimate_normals() first");
        return normals[n];
    }

    Point get_triangle_normal(int n) const {
        if (n >= (int) triangles.size() || (int) triangles.size() == 0)
            throw MeshException("get_triangle: index exceeds face dimensions");
        return triangles[n].normal();
    }

    int get_triangle_vertexID(int n, int i) const { return triangles[n].get_vertex_no(i); }
    std::vector<std::shared_ptr<Mpoint>> get_points() const { return points; }

    int get_total_neighbours(int i) const {
        return points.at(i)->nneighbours();
    }

    std::vector<std::vector<double>> get_face_angles() const {

        std::vector<std::vector<double>> all_angles(triangles.size());
        for (const auto& triangle : triangles)
            all_angles.emplace_back(triangle.get_angles());

        return all_angles;
    }

    Point local_normal(int pt) const;
    void estimate_normals();
    int get_resolution() const;
    double calculate_MeanVD() const;
    double calculate_MaxVD() const;

    //---MANIPULATION---//
    std::vector<std::shared_ptr<Mpoint>>& get_all_points() { return points; }
    inline void set_coord(int i, const Point &p) { points[i]->set_coord(p); }
    void set_pvalue(unsigned int i, float val, int dim = 0);
    void set_pvalues(const NEWMAT::Matrix &M, bool appendFieldData = false);
    void initialize_pvalues(int dim = 1, bool appendFieldData = false); // initialize new data array
    void clear_triangles() { triangles.clear(); }
    void push_point(const std::shared_ptr<Mpoint>& mesh_point) { points.push_back(mesh_point); }
    void push_pvalues(const std::vector<float>& pvals) { pvalues.push_back(pvals); }
    void push_tvalues(const std::vector<float>& tvals) { tvalues.push_back(tvals); }
    void push_triangle(const Triangle &t);
    void clear();
    void clear_data();
    Point estimate_origin() const;

    //---ITERATORS---//
    std::vector<std::shared_ptr<Mpoint>>::const_iterator vbegin() const { return points.begin(); };
    std::vector<std::shared_ptr<Mpoint>>::const_iterator vend() const { return points.end(); };
    std::vector<Triangle>::iterator tbegin() { return triangles.begin(); };
    std::vector<Triangle>::const_iterator tbegin() const { return triangles.begin(); };
    std::vector<Triangle>::iterator tend() { return triangles.end(); };
    std::vector<Triangle>::const_iterator tend() const { return triangles.end(); };
    std::vector<int>::const_iterator nbegin(int i) const { return points[i]->get_nID().begin(); };
    std::vector<int>::const_iterator nend(int i) const { return points[i]->get_nID().end(); };
    std::vector<int>::const_iterator tIDbegin(int i) const { return points[i]->get_trID().begin(); };
    std::vector<int>::const_iterator tIDend(int i) const { return points[i]->get_trID().end(); };

    //---LOAD/WRITE FUNCTIONS---//
    enum class FileType { ASCII, VTK, GIFTI, MATRIX, DPV, DEFAULT };

    Mesh::FileType meshFileType(const std::string &filename) const;

    void load(const std::string &filename, bool loadSurfaceData = true, bool appendFieldData = false);
    void load_gifti(const std::string &filename, bool loadSurfaceData = true, bool appendFieldData = false);
    void load_ascii(const std::string &filename, bool loadSurfaceData = true, bool appendFieldData = false);
    void load_vtk(const std::string &filename);
    void load_matrix(const std::string &filename, const Mesh::FileType &type);
    void load_ascii_file(const std::string &filename);

    void save(const std::string &filename) const;
    void save_ascii(const std::string &s) const;
    void save_gifti(const std::string &s) const;
    void save_vtk(const std::string &s) const;
    void save_dpv(const std::string &s) const;
    void save_matrix(const std::string &s) const;

    std::vector<float> getPointsAsVectors() const;
    std::vector<float> getValuesAsVectors(int d = 0) const;
    std::vector<int> getTrianglesAsVector() const;
};

//---Functions to generate regular icosahedron---//
Mesh make_mesh_from_icosa(int n);
void retessellate(Mesh&);
void retessellate(Mesh&, std::vector<std::vector<int>>&);

//---Utilities---//
void check_scale(Mesh& in, const Mesh& ref);
void true_rescale(Mesh& m, double rad);
void recentre(Mesh& sphere);
Mesh create_exclusion(const Mesh& data_mesh, float thrl, float thru);
double compute_vertex_area(int, const Mesh &); // averages adjoining face areas for each vertex

//---Mesh operators---//
bool operator==(const Mesh &M1, const Mesh &M2);

} //namespace newresampler

#endif  //NEWRESAMPLER_MESH_H
