# gMSM - groupwise registration mode of MSM

These scripts demonstrate the usage of gMSM. In order to run gMSM set the followings.

## Usage

Create a directory named input and copy your cohort data files and mesh files into a directory with your cohort ID. Also, copy the config file and template surface file.
Set the directory structure as shown in the following example using the tree command line tool:
```console
$ tree .
.
└── input
    └── 2133
        ├── data
        │   ├── 117122.L.sulc.affine.ico6.shape.gii
        │   ├── 117123.L.sulc.affine.ico6.shape.gii
        │   ├── 128026.L.sulc.affine.ico6.shape.gii
        │   ├── 139637.L.sulc.affine.ico6.shape.gii
        │   ├── 561444.L.sulc.affine.ico6.shape.gii
        │   ├── 660951.L.sulc.affine.ico6.shape.gii
        │   └── 677968.L.sulc.affine.ico6.shape.gii
        ├── gMSM_config.txt
        ├── meshes
        │   ├── 117122.surf.gii
        │   ├── 117123.surf.gii
        │   ├── 128026.surf.gii
        │   ├── 139637.surf.gii
        │   ├── 561444.surf.gii
        │   ├── 660951.surf.gii
        │   └── 677968.surf.gii
        └── sunet.ico-6.template.surf.gii
```

Set the first few lines of run_gMSM.sh and stats.py according to your settings and cohort.
