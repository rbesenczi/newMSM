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
#include "mpoint.h"

namespace newresampler {

bool Mpoint::is_triangle(int j) {
    auto it = std::find(trID.begin(), trID.end(), j);
    return it != trID.end();
}

bool Mpoint::is_neighbour(int j) {
    auto it = std::find(nID.begin(), nID.end(), j);
    return it != nID.end();
}

bool operator==(const Mpoint &p1, const Mpoint &p2) {
    return p1.get_coord() == p2.get_coord();
}

bool operator!=(const Mpoint &p1, const Mpoint &p2) {
    return !(p1==p2);
}

bool operator==(const Mpoint &p1, const Point &p2) {
    return p1.get_coord() == p2;
}

Point operator-(const Mpoint &p1, const Mpoint &p2) {
    return p1.get_coord() - p2.get_coord();
}

Point operator-(const Point &p1, const Mpoint &p2) {
    return p1 - p2.get_coord();
}

} // namespace newresampler
