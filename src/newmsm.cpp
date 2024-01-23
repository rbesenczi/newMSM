#include "msmOptions.h"
#include "NewMeshReg/group_mesh_registration.h"
#include "NewMeshReg/group_coregistration.h"

msmOptions *msmOptions::gopt = nullptr;

int main(int argc, char *argv[]) try {

    std::cout << "This is newMSM v0.5.1-BETA." << std::endl;

    msmOptions &opts = msmOptions::getInstance();
    opts.parse_command_line(argc, argv);

    if(opts.in_register.set())
        std::cout << "Warning! --in_register parameter is removed from newMSM." << std::endl;

    if(opts.multiresolutionlevels.set())
        std::cout << "Warning! --levels parameter is removed from newMSM." << std::endl;

    if(opts.smoothoutput.set())
        std::cout << "Warning! --smoothout parameter is removed from newMSM." << std::endl;

    if(opts.groupwise.value())
    {
        newmeshreg::Group_Mesh_registration GMR;
        if (opts.printoptions.value()) GMR.print_config_options();

        if (opts.verbose.value()) GMR.set_verbosity(opts.verbose.value());
        if (opts.debug.value()) GMR.set_debug(opts.debug.value());
        GMR.set_outdir(opts.outbase.value());

        GMR.set_inputs(opts.meshes.value());
        GMR.set_template(opts.templatemesh.value());
        GMR.set_data_list(opts.data.value());

        GMR.run_multiresolutions(opts.parameters.value());
    }
    if(opts.cogroup.value())
    {
        std::cout << "coGroup mode selected" << std::endl;
        newmeshreg::Group_coregistration GCR;
        if (opts.printoptions.value()) GCR.print_config_options();
        if (opts.verbose.value()) GCR.set_verbosity(opts.verbose.value());
        if (opts.debug.value()) GCR.set_debug(opts.debug.value());
        GCR.set_outdir(opts.outbase.value());

        GCR.set_warps(opts.meshesA.value(), opts.meshesB.value());
        GCR.set_meshes(opts.meshA.value(), opts.meshB.value());

        GCR.set_data(opts.meanA.value(), opts.meanB.value());
        GCR.set_template(opts.templatemesh.value());

        //GCR.run_multiresolutions(opts.parameters.value());
    }
    else
    {
        newmeshreg::Mesh_registration MR;
        if (opts.printoptions.value()) MR.print_config_options();

        if (opts.verbose.value()) MR.set_verbosity(opts.verbose.value());
        if (opts.debug.value()) MR.set_debug(opts.debug.value());
        MR.set_outdir(opts.outbase.value());

        MR.set_input(opts.inputmesh.value());
        if (opts.referencemesh.value().empty()) MR.set_reference(opts.inputmesh.value());
        else MR.set_reference(opts.referencemesh.value());

        if (!opts.inputanatmesh.value().empty()) {
            if (opts.referenceanatmesh.value().empty())
                throw newmeshreg::MeshregException("Error: must supply both anatomical meshes or none");

            MR.set_anatomical(opts.inputanatmesh.value(), opts.referenceanatmesh.value());
        }

        MR.set_output_format(opts.outformat.value());

        if (opts.transformed_sphere.set()) MR.set_transformed(opts.transformed_sphere.value());
        if (opts.cfweight_in.set()) MR.set_input_cfweighting(opts.cfweight_in.value());
        if (opts.cfweight_ref.set()) MR.set_reference_cfweighting(opts.cfweight_ref.value());

        MR.set_CMpathin(opts.CMmatrixin.value());
        MR.set_CMpathref(opts.CMmatrixref.value());

        MR.run_multiresolutions(opts.parameters.value());
    }

    return EXIT_SUCCESS;

} catch (newmeshreg::MeshregException& e) {
    e.what();
    exit(EXIT_FAILURE);
} catch (std::exception& e) {
    e.what();
    exit(EXIT_FAILURE);
}
