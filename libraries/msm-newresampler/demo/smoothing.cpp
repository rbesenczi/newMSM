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
#include "miscmaths/bfmatrix.h"
#include <boost/program_options.hpp>
#include <memory>

#define RAD 100

void smooth(const std::string& in_sphere_name, const std::string& data_name, 
        double sigma, const std::string& output_name) {
    //tests data smoothing on the same mesh

    newresampler::Mesh in_sphere, smoothed;
    in_sphere.load(in_sphere_name, true, false);
    in_sphere.load(data_name, false, false);

    newresampler::true_rescale(in_sphere, RAD);

    smoothed = newresampler::smooth_data(in_sphere, in_sphere, sigma);

    smoothed.save(output_name + "-smoothed_data.func");
}

int main(int argc, char **argv)
try {

    boost::program_options::options_description desc("Arguments");
    desc.add_options()
            ("help", "This message")
            ("metric_in", boost::program_options::value<std::string>(), "metric file to resample")
            ("current_sphere", boost::program_options::value<std::string>(), "the current sphere the metric file is on")
            ("sigma", boost::program_options::value<double>(), "sigma parameter of smoothing (e.g. 0.5)")
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

    if (vm.count("sigma")) {
#ifdef DEBUG
        std::cout << "sigma was set to "
                  << vm["sigma"].as<double>() << ".\n";
#endif
    } else {
        std::cout << "sigma was not set, but required.\n";
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

    smooth(vm["current_sphere"].as<std::string>(), vm["metric_in"].as<std::string>(), vm["sigma"].as<double>(), vm["output"].as<std::string>());

    return 0;

} catch (std::exception& some_exception) {
    std::cout << some_exception.what() << '\n';
    return 1;
}
