# msm-newmeshreg library for MSM software

This is a rework and extension of MSM mesh registration library. The code and licenses of the original version can be found at [this link](https://github.com/ecr05/MSM_HOCR). The user guide of that version can be found [here](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/MSM).

## Major updates

Most of the code from the original version has been replaced or revised, but some were reused. Most of the updates affect performance:
 - application of octree search during registration
 - parallelisation of cost calculation during registration

Please note, the current version is 0.5.1-BETA. All feedbacks are much appreciated.

## Install from source

1. Download and install FSL. For more information, please visit the FSL [webpage](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/).

2. Download and install sources and demo applications.
    ```console
    git clone https://github.com/rbesenczi/newMSM.git
    cd newMSM/libraries/msm-newmeshreg/src/
    make && make install
    ```

Please note the /licenses folder - including copy of the ELC and FastPD licence info. (for information only; I supply the ELC library with permission from the author (for research use only)).