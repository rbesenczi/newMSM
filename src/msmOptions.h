#ifndef MSMOPTIONS_H
#define MSMOPTIONS_H

#include "utils/options.h"

class msmOptions {

public:
    static msmOptions &getInstance();

    ~msmOptions() { delete gopt; }

    Utilities::Option<bool> help;
    Utilities::Option<bool> verbose;
    Utilities::Option<bool> printoptions;
    Utilities::Option<bool> debug;
    Utilities::Option<std::string> inputmesh;
    Utilities::Option<std::string> referencemesh;
    Utilities::Option<std::string> inputanatmesh;
    Utilities::Option<std::string> referenceanatmesh;
    Utilities::Option<std::string> CMmatrixin;
    Utilities::Option<std::string> CMmatrixref;
    Utilities::Option<std::string> transformed_sphere;
    Utilities::Option<std::string> cfweight_in;
    Utilities::Option<std::string> cfweight_ref;
    Utilities::Option<std::string> outbase;
    Utilities::Option<std::string> outformat;
    Utilities::Option<std::string> parameters;

    //Groupwise params
    Utilities::Option<bool> groupwise;
    Utilities::Option<std::string> meshes;
    Utilities::Option<std::string> templatemesh;
    Utilities::Option<std::string> data;

    //Co-groupwise params
    Utilities::Option<bool> cogroup;
    Utilities::Option<std::string> meanA;
    Utilities::Option<std::string> meanB;
    Utilities::Option<std::string> meshA;
    Utilities::Option<std::string> meshB;
    Utilities::Option<std::string> meshesA;
    Utilities::Option<std::string> meshesB;

    // Removed params
    Utilities::Option<std::string> in_register;
    Utilities::Option<int> multiresolutionlevels;
    Utilities::Option<float> smoothoutput;

    bool parse_command_line(int argc, char **argv);

    const msmOptions& operator=(const msmOptions&) = delete;
    msmOptions(const msmOptions&) = delete;
    const msmOptions& operator=(msmOptions&&) = delete;
    msmOptions(msmOptions&&) = delete;

private:
    msmOptions();   //singleton class

    Utilities::OptionParser options;
    static msmOptions* gopt;
};

inline msmOptions& msmOptions::getInstance() {

    if (gopt == nullptr)
        gopt = new msmOptions();

    return* gopt;
}

inline msmOptions::msmOptions() :
        options("newmsm", "newmsm [options]\n"),
        help(std::string("-h,--help"), false,
             std::string("display this message"),
             false, Utilities::no_argument),
        verbose(std::string("-v,--verbose"), false,
                std::string("switch on diagnostic messages"),
                false, Utilities::no_argument),
        printoptions(std::string("-p,--printoptions"), false,
                     std::string("print configuration file options"),
                     false, Utilities::no_argument),
        debug(std::string("--debug"), false,
              std::string("run debugging or optimising options"),
              false, Utilities::no_argument,false),
        groupwise(std::string("--groupwise"), false,
              std::string("Run newMSM in groupwise mode"),
              false, Utilities::no_argument,false),
        meshes(std::string("--meshes"), std::string(""),
               std::string("list of paths to input meshes (available formats: VTK, ASCII, GIFTI). Needs to be a sphere"),
               false , Utilities::requires_argument),
        templatemesh(std::string("--template"), std::string(""),
                     std::string("templates sphere for resampling (available formats: VTK, ASCII, GIFTI). Needs to be a sphere"),
                     false , Utilities::requires_argument),
        cogroup(std::string("--cogroup"), false,
                  std::string("Run newMSM in group co-registering mode"),
                  false, Utilities::no_argument,false),
        data(std::string("--data"), std::string(""),
             std::string("list of paths to the data"),
             false , Utilities::requires_argument),
        meanA(std::string("--meanA"), std::string(""),
              std::string("Mean of group A"),
              false , Utilities::requires_argument),
        meanB(std::string("--meanB"), std::string(""),
              std::string("Mean of group B"),
              false , Utilities::requires_argument),
        meshA(std::string("--meshA"), std::string(""),
              std::string("Mesh of mean of group A"),
              false , Utilities::requires_argument),
        meshB(std::string("--meshB"), std::string(""),
              std::string("Mesh of mean of group B"),
              false , Utilities::requires_argument),
        meshesA(std::string("--meshesA"), std::string(""),
              std::string("List of warps of group A"),
              false , Utilities::requires_argument),
        meshesB(std::string("--meshesB"), std::string(""),
              std::string("List of warps of group B"),
              false , Utilities::requires_argument),
        inputmesh(std::string("--inmesh"), std::string(""),
                  std::string("input mesh (available formats: VTK, ASCII, GIFTI). Needs to be a sphere"),
                  false , Utilities::requires_argument),
        referencemesh(std::string("--refmesh"), std::string(""),
                      std::string("reference mesh (available formats: VTK, ASCII, GIFTI). Needs to be a sphere. If not included algorithm assumes reference mesh is equivalent input"),
                      false , Utilities::requires_argument),
        inputanatmesh(std::string("--inanat"), std::string(""),
                      std::string("input anatomical mesh (available formats: VTK, ASCII, GIFTI). For Example, white, pial, midthickness  (must either supply both input and reference anatomical surfaces or none)"),
                      false , Utilities::requires_argument,false),
        referenceanatmesh(std::string("--refanat"), std::string(""),
                          std::string("reference mesh (available formats: VTK, ASCII, GIFTI). For Example, white, pial, midthickness (must either supply both input and reference anatomical surfaces or none)"),
                          false , Utilities::requires_argument,false),
        CMmatrixin(std::string("--indata"), std::string(""),
                   std::string("scalar or multivariate data for input - can be ASCII (.asc,.dpv,.txt) or GIFTI (.func.gii or .shape.gii) "),
                   false, Utilities::requires_argument),
        CMmatrixref(std::string("--refdata"), std::string(""),
                    std::string("scalar or multivariate data for reference - can be ASCII (.asc,.dpv,.txt) or GIFTI (.func.gii or .shape.gii) "),
                    false, Utilities::requires_argument),
        transformed_sphere(std::string("--trans"), std::string(""),
                           std::string("\t Transformed source mesh (output of a previous registration). Use this to initiliase the current registration."),
                           false , Utilities::requires_argument),
        cfweight_in(std::string("--inweight"), std::string(""),
                    std::string("cost function weighting for input - weights data in these vertices when calculating similarity (ASCII or GIFTI). Can be multivariate provided dimension equals that of data "),
                    false , Utilities::requires_argument),
        cfweight_ref(std::string("--refweight"), std::string(""),
                     std::string("cost function weighting for reference - weights data in these vertices when calculating similarity (ASCII or GIFTI). Can be multivariate provided dimension equals that of data"),
                     false , Utilities::requires_argument),
        outbase(std::string("-o,--out"), std::string(""),
                std::string("output basename"),
                false, Utilities::requires_argument),
        in_register(std::string("--in_register"), std::string(""),
                   std::string("\t Warning! This option is removed from newMSM."),
                   false, Utilities::requires_argument),
        multiresolutionlevels(std::string("--levels"),0,
                    std::string("\t Warning! This option is removed from newMSM."),
                    false, Utilities::requires_argument),
        smoothoutput(std::string("--smoothout"), 0.0,
                   std::string("\t Warning! This option is removed from newMSM."),
                   false, Utilities::requires_argument),
        outformat(std::string("-f,--format"), std::string("GIFTI"),
                  std::string("format of output files, can be: GIFTI, VTK, ASCII or ASCII_MAT (for full details of output file formats see MSM wiki)"),
                  false, Utilities::requires_argument),
        parameters(std::string("--conf"), std::string(""),
                   std::string("\tconfiguration file "),
                   false, Utilities::requires_argument)
    {
    try {
        options.add(help);
        options.add(verbose);
        options.add(printoptions);
        options.add(debug);
        options.add(groupwise);
        options.add(meshes);
        options.add(templatemesh);
        options.add(data);
        options.add(cogroup);
        options.add(meanA);
        options.add(meanB);
        options.add(meshA);
        options.add(meshB);
        options.add(meshesA);
        options.add(meshesB);
        options.add(inputmesh);
        options.add(referencemesh);
        options.add(inputanatmesh);
        options.add(referenceanatmesh);
        options.add(CMmatrixin);
        options.add(CMmatrixref);
        options.add(transformed_sphere);
        options.add(cfweight_in);
        options.add(cfweight_ref);
        options.add(outbase);
        options.add(outformat);
        options.add(parameters);
        options.add(in_register);
        options.add(multiresolutionlevels);
        options.add(smoothoutput);
    }
    catch(Utilities::X_OptionError& e)
    {
        options.usage();
        std::cerr << '\n' << e.what() << std::endl;
    }
    catch(std::exception &e)
    {
        std::cerr << e.what() << std::endl;
    }
}

inline bool msmOptions::parse_command_line(int argc, char **argv) {

    for (unsigned int a = options.parse_command_line(argc, argv); a < argc; a++);

    if (!(printoptions.value()))
    {
        if (help.value())
        {
            options.usage();
            exit(EXIT_SUCCESS);
        }

        if (!options.check_compulsory_arguments())
        {
            options.usage();
            exit(EXIT_FAILURE);
        }
    }

    return true;
}

#endif // MSMOPTIONS_H
