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
#include "newresampler/resampler.h"
#include <boost/program_options.hpp>

#define RAD 100

void resample_mesh(const std::string& in_sphere_name, const std::string& data_name, 
        int ico_dim, const std::string& output_name) {

    newresampler::Mesh in_sphere;
    in_sphere.load(in_sphere_name, true, false);
    in_sphere.load(data_name, false, false);

    newresampler::Mesh ico = newresampler::make_mesh_from_icosa(ico_dim);

    newresampler::true_rescale(in_sphere,RAD);
    newresampler::true_rescale(ico, RAD);

    newresampler::Mesh resampled = newresampler::metric_resample(in_sphere, ico);
    
    resampled.save(output_name + "-resampled_data.func");
}

int main(int argc, char **argv)
try {

    boost::program_options::options_description desc("Arguments");
    desc.add_options()
            ("help", "This message")
            ("metric_in", boost::program_options::value<std::string>(), "metric file to resample")
            ("current_sphere", boost::program_options::value<std::string>(), "the current sphere the metric file is on")
            ("ico", boost::program_options::value<int>(), "the order of the generated regular icosahedron to resample to (3-6)")
            ("output", boost::program_options::value<std::string>(), "output base file name")
            ;

    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    boost::program_options::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 0;
    }

    if (vm.count("metric_in")) {
#ifdef DEBUG
        std::cout << "metric_in was set to "
                  << vm["data"].as<std::string>() << ".\n";
#endif
    } else {
        std::cout << "metric_in was not set, but required.\n";
        return 1;
    }

    if (vm.count("current_sphere")) {
#ifdef DEBUG
        std::cout << "current_sphere was set to " << vm["current_sphere"].as<std::string>() << ".\n";
#endif
    } else {
        std::cout << "current_sphere was not set, but required.\n";
        std::cout << desc << '\n';
        return 1;
    }

    if (vm.count("ico")) {
#ifdef DEBUG
        std::cout << "ico was set to "
                  << vm["ico"].as<int>() << ".\n";
#endif
        if (vm["ico"].as<int>() > 6 || vm["ico"].as<int>() < 2)
        {
            std::cout << desc << '\n';
            throw std::runtime_error("Invalid ico dimension");
        }
    } else {
        std::cout << "ico was not set, but required.\n";
        std::cout << desc << '\n';
        return 1;
    }

    if (vm.count("output")) {
#ifdef DEBUG
        std::cout << "output was set to "
                  << vm["output"].as<std::string>() << ".\n";
#endif
    } else {
        std::cout << "output was not set, but required.\n";
        std::cout << desc << '\n';
        return 1;
    }

    resample_mesh(vm["current_sphere"].as<std::string>(), vm["metric_in"].as<std::string>(),
        vm["ico"].as<int>(), vm["output"].as<std::string>());

    return 0;

} catch (std::exception& some_exception) {
    std::cout << some_exception.what() << '\n';
    return 1;
}
