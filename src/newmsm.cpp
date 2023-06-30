#include "msmOptions.h"
#include "NewMeshReg/mesh_registration.h"

msmOptions *msmOptions::gopt = nullptr;

int main(int argc, char *argv[]) {

    std::cout << "This is newMSM v0.4.1-BETA." << std::endl;

    msmOptions &opts = msmOptions::getInstance();
    opts.parse_command_line(argc, argv);

    newmeshreg::Mesh_registration MR;

    if(opts.in_register.set())
        std::cout << "Warning! --in_register parameter is removed from newMSM." << std::endl;

    if(opts.multiresolutionlevels.set())
        std::cout << "Warning! --levels parameter is removed from newMSM." << std::endl;

    if(opts.smoothoutput.set())
        std::cout << "Warning! --smoothout parameter is removed from newMSM." << std::endl;

    if (opts.printoptions.value())
    {
        MR.print_config_options();
        return 0;
    }

    if (opts.verbose.value()) MR.set_verbosity(opts.verbose.value());
    if (opts.debug.value()) MR.set_debug(opts.debug.value());

    MR.set_input(opts.inputmesh.value());

    if (opts.referencemesh.value().empty()) MR.set_reference(opts.inputmesh.value());
    else MR.set_reference(opts.referencemesh.value());

    if (!opts.inputanatmesh.value().empty())
    {
        if (opts.referenceanatmesh.value().empty())
        {
            std::cout << " Error: must supply both anatomical meshes or none " << std::endl;
            exit(1);
        }
        MR.set_anatomical(opts.inputanatmesh.value(), opts.referenceanatmesh.value());
    }

    MR.set_outdir(opts.outbase.value());
    MR.set_output_format(opts.outformat.value());

    if (opts.transformed_sphere.set()) MR.set_transformed(opts.transformed_sphere.value());
    if (opts.cfweight_in.set()) MR.set_input_cfweighting(opts.cfweight_in.value());
    if (opts.cfweight_ref.set()) MR.set_reference_cfweighting(opts.cfweight_ref.value());

    MR.set_CMpathin(opts.CMmatrixin.value());
    MR.set_CMpathref(opts.CMmatrixref.value());

    MR.run_multiresolutions(opts.parameters.value());

    return 0;
}
