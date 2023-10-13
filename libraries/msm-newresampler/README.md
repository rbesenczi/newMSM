# msm-newresampler library for MSM software

This is a rework and extension of MSM resampler and mesh library. The original library can be found at [this link](https://git.fmrib.ox.ac.uk/fsl/newmesh/-/tree/master).

## Major updates

Most of the code from the original version has been replaced or revised, but some were reused. Most of the updates affect performance:

- complete rework of the mesh resampling functionality with an octree based nearest triangle search algorithm (inspired by the wb_command of Connectome Workbench)
- parallelisation of mesh resampling
- reorganisation of source code and standard code cleanup/code revision.

Please note, the current version is 0.4.2-BETA. All feedbacks are much appreciated.

A few demo application can be found in the demo folder. Please see installation and usage instructions below.

## Install from source

1. Download and install FSL. For more information, please visit the FSL [webpage](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/).

2. Download and install sources and demo applications.
    ```console
    git clone https://github.com/rbesenczi/newMSM.git
    cd newMSM/libraries/msm-newresampler/src/
    git checkout tags/v0.4.2
    make && make install
    cd ../demo
    make && make install
    ```
