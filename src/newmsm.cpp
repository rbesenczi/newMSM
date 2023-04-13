#include "msmOptions.h"
#include "NewMeshReg/mesh_registration.h"

msmOptions *msmOptions::gopt = nullptr;

int main(int argc, char *argv[]) {

    std::cout << "This is newMSM v.0.1.1-BETA." << std::endl;

    msmOptions &opts = msmOptions::getInstance();
    Utilities::Log &logger = Utilities::LogSingleton::getInstance();
    std::string CMpathin;

    opts.parse_command_line(argc, argv, logger);
    newmeshreg::Mesh_registration MR;

    //---INITIALISE VARIABLES---//
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

    // add anatomical meshes for stress and strain
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

    if (opts.in_register.set())
    {
        // if data is supplied at a different resolution to the input or
        // transformed mesh it is necessary to resample i.e. HCP 32K to native
        // the in_register sphere HAS to be the sphere on which the data
        // is supplied and the input or transformed mesh MUST be in alignment with it.
        newresampler::Mesh in_register, target;

        if (opts.transformed_sphere.set()) target.load(opts.transformed_sphere.value());
        else target.load(opts.inputmesh.value());
        in_register.load(opts.in_register.value());

        std::shared_ptr<MISCMATHS::BFMatrix> DATA;
        newmeshreg::set_data(opts.CMmatrixin.value(), DATA, in_register);

        in_register.set_pvalues(DATA->AsMatrix());
        in_register = newresampler::metric_resample(in_register, target);
        target.set_pvalues(DATA->AsMatrix());

        std::string fname = opts.outbase.value() + "input_in_register.func.gii";
        target.save(fname);
        CMpathin = fname;
    }
    else
        CMpathin = opts.CMmatrixin.value();

    MR.set_CMpathin(CMpathin);
    MR.set_CMpathref(opts.CMmatrixref.value());

    MR.run_multiresolutions(opts.multiresolutionlevels.value(), opts.smoothoutput.value(), opts.parameters.value());

    if (opts.in_register.set())
    {
        // if data is supplied at a lower mesh resolution,
        // then resample final warp accordingly.
        // This is typically a HCP formating issue.
        newresampler::Mesh ORIG, FINAL, in_register;

        if (opts.transformed_sphere.set()) ORIG.load(opts.transformed_sphere.value());
        else ORIG.load(opts.inputmesh.value());

        in_register.load(opts.in_register.value());
        FINAL.load(opts.outbase.value() + "sphere.reg" + MR.get_surf_format());
        in_register = newresampler::surface_resample(in_register, ORIG, FINAL);
        in_register.save(opts.outbase.value() + "sphere.in_register.reg." + MR.get_surf_format());
    }

    return 0;
}
