# newMSM - the new version of the Multimodal Surface Matching software

## Download and install
NewMSM is available to download from FSL or Github. Please follow the steps described below. (Note: tested on Ubuntu 20.04, 22.04, 24.04 and MacOS Sonoma 14.7.1.)

### Download and install binary

1. Download and install FSL. For more information, please visit the FSL [webpage](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/).
2. Download and install newMSM from the FSL developer branch by executing the following command on Terminal:
```console
fsl/bin/conda create --override-channels -c https://fsl.fmrib.ox.ac.uk/fsldownloads/fslconda/development/ -c https://fsl.fmrib.ox.ac.uk/fsldownloads/fslconda/public/ -c conda-forge -p ./msm-env fsl-newmsm
```
3. Activate this environment:
```console
conda activate ./msm-env/
```
4. You can test with the `newmsm -p` command.

### Install from source

1. Download and install FSL. For more information, please visit the FSL [webpage](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/).
2. Create a directoy where you can download sources:
    ```console
    mkdir newmsm_sources
    cd newmsm_sources
    ```
3. Download and install newMSM sources.
    1. newMSM application and dependencies
        ```console
        git clone https://github.com/rbesenczi/newMSM.git
        cd newMSM/libraries/msm-newresampler/src/
        make && make install
        cd ../../msm-newmeshreg/src/
        make && make install
        cd ../../../src/
        make && make install
        ```
3. The binary `newmsm` should be in the folder `${FSLDEVDIR}/bin`.
4. You can test it with the `${FSLDEVDIR}/bin/newmsm -p` command.

Most features remained the same as it was in the original MSM implementation. A detailed user guide can be found [here](https://github.com/rbesenczi/newMSM/docs/guide.md).

## Links

FSL download and install [instructions](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation).
