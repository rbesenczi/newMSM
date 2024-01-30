/*
Copyright (c) 2024 King's College London, MeTrICS Lab, Renato Besenczi

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
#include "newresampler/resampler.h"
#include <boost/program_options.hpp>

void applywarp(const std::string& to_be_def, const std::string& orig_sph,
               const std::string& warp_sph, const std::string& output_name) {

    newresampler::Mesh to_be_deformed, original_sphere, warp;
    to_be_deformed.load(to_be_def, true, false);
    original_sphere.load(orig_sph, true, false);
    warp.load(warp_sph, true, false);

    newresampler::true_rescale(to_be_deformed,RAD);
    newresampler::true_rescale(original_sphere, RAD);
    newresampler::true_rescale(warp, RAD);

    newresampler::barycentric_mesh_interpolation(to_be_deformed, original_sphere, warp);

    to_be_deformed.save(output_name + "-warped.surf");
}

int main(int argc, char **argv)
try {

    boost::program_options::options_description desc("Arguments");
    desc.add_options()
            ("help", "This message")
            ("to_be_deformed", boost::program_options::value<std::string>(), "mesh to be deformed")
            ("original_sphere", boost::program_options::value<std::string>(), "the original sphere the mesh is on")
            ("warp", boost::program_options::value<std::string>(), "the warp to be applied")
            ("output", boost::program_options::value<std::string>(), "output base file name")
            ;

    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    boost::program_options::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        exit(EXIT_SUCCESS);
    }

    if (!vm.count("to_be_deformed")) {
        std::cout << "to_be_deformed was not set, but required.\n";
        exit(EXIT_FAILURE);
    }

    if (!vm.count("original_sphere")) {
        std::cout << "original_sphere was not set, but required.\n";
        std::cout << desc << '\n';
        exit(EXIT_FAILURE);
    }

    if (!vm.count("warp")) {
        std::cout << "warp was not set, but required.\n";
        std::cout << desc << '\n';
        exit(EXIT_FAILURE);
    }

    if (!vm.count("output")) {
        std::cout << "output was not set, but required.\n";
        std::cout << desc << '\n';
        exit(EXIT_FAILURE);
    }

    applywarp(vm["to_be_deformed"].as<std::string>(), vm["original_sphere"].as<std::string>(),
              vm["warp"].as<std::string>(), vm["output"].as<std::string>());

    exit(EXIT_SUCCESS);

} catch (std::exception& some_exception) {
    std::cout << some_exception.what() << '\n';
    exit(EXIT_FAILURE);
}
