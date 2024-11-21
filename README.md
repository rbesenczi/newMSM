# newMSM - the new version of the Multimodal Surface Matching software

Multimodal Surface Matching (MSM) is an advanced tool for cortical surface mapping, which has been shown to outperform competing methods for the alignment of multimodal data. MSM has been developed and tested using FreeSurfer extracted surfaces, however, in principle the tool will work with any cortical surface extraction method provided the surfaces can be mapped to the sphere. The key advantage of the method is that alignment may be driven using a wide variety of univariate (sulcal depth, curvature, myelin), multivariate (Task fMRI, or Resting State Networks) or multimodal (combinations of folding, myelin and fMRI) feature sets.

The main MSM tool is currently run from the command line using the program `newmsm`. This enables fast alignment of spherical cortical surfaces by utilising discrete optimisation framework, which significantly reduces the search space of possible deformations for each vertex, and allows flexibility with regards to the choice of similarity metric used to match the images.

NewMSM is a new implementation with several improvements that made the MSM method computationally more efficient. Improvements include a completely re-implemented mesh resampling library with a new nearest neighbour search algorithm (octree search), an option for multicore CPU utilisation and several others. In general, the average runtime of a registration process is now 25% of that of the original MSM implementation (and 5% using multicore CPU utilisation). The operation of newMSM is supposed to be the same as the previous MSM implementation. We notify the user about any changes that have been made in adequate command line messages. NewMSM now contains an implementation of a groupwise surface registration method that is described later in this guide.

Before use, read our [user's guide](docs/guide.md).

Installation instruction can be found [here](docs/install.md).

The original MSM implementation can be found [here](https://github.com/ecr05/MSM_HOCR/).

## Major updates

Most of the code from the original version has been replaced, but some were reused after a major revision. Most of the updates affect performance:
 - complete rework of the mesh resampling functionality with an octree based nearest triangle search algorithm (inspired by the wb_command of Connectome Workbench)
 - parallelisation of mesh resampling
 - application of octree search during registration
 - parallelisation of cost calculation during registration

As a new feature, groupwise cortical surface registration is now available. Please see our [user's guide](https://github.com/rbesenczi/newMSM/docs/guide.md) for more information.

Please note, the current version is 0.8.0-BETA. All feedbacks are much appreciated.

## Referencing

If you wish to use MSM, please reference our papers in any resulting publication.

> Emma C. Robinson, Saad Jbabdi, Matthew F. Glasser, Jesper Andersson, Gregory C. Burgess, Michael P. Harms, Stephen M. Smith, David C. Van Essen, Mark Jenkinson, MSM: A new flexible framework for Multimodal Surface Matching, NeuroImage, Volume 100, 2014, Pages 414-426, ISSN 1053-8119, [https://doi.org/10.1016/j.neuroimage.2014.05.069](https://doi.org/10.1016/j.neuroimage.2014.05.069).


> N. Komodakis and G. Tziritas, Approximate Labeling via Graph Cuts Based on Linear Programming, in IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 29, no. 8, pp. 1436-1453, Aug. 2007, [https://doi.org/10.1109/TPAMI.2007.1061](https://doi.org/10.1109/TPAMI.2007.1061). 

For those using the HCP pipelines and/or MSM with Higher Order Smoothness Constraints, please also reference:

> Emma C. Robinson, Kara Garcia, Matthew F. Glasser, Zhengdao Chen, Timothy S. Coalson, Antonios Makropoulos, Jelena Bozek, Robert Wright, Andreas Schuh, Matthew Webster, Jana Hutter, Anthony Price, Lucilio Cordero Grande, Emer Hughes, Nora Tusor, Philip V. Bayly, David C. Van Essen, Stephen M. Smith, A. David Edwards, Joseph Hajnal, Mark Jenkinson, Ben Glocker, Daniel Rueckert, Multimodal surface matching with higher-order smoothness constraints, NeuroImage, Volume 167, 2018, Pages 453-465, ISSN 1053-8119, [https://doi.org/10.1016/j.neuroimage.2017.10.037](https://doi.org/10.1016/j.neuroimage.2017.10.037).

> Ishikawa, Hiroshi. Higher-order clique reduction without auxiliary variables, In Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition, pp. 1362-1369. 2014.

> Matthew F. Glasser, Stamatios N. Sotiropoulos, J. Anthony Wilson, Timothy S. Coalson, Bruce Fischl, Jesper L. Andersson, Junqian Xu, Saad Jbabdi, Matthew Webster, Jonathan R. Polimeni, David C. Van Essen, Mark Jenkinson, The minimal preprocessing pipelines for the Human Connectome Project, NeuroImage, Volume 80, 2013, Pages 105-124, ISSN 1053-8119, [https://doi.org/10.1016/j.neuroimage.2013.04.127](https://doi.org/10.1016/j.neuroimage.2013.04.127).

If you use the groupwise method, please also reference:

> Besenczi, R., Guo, Y., Robinson, E.C. (2024). High Performance Groupwise Cortical Surface Registration with Multimodal Surface Matching. In: Modat, M., Simpson, I., Špiclin, Ž., Bastiaansen, W., Hering, A., Mok, T.C.W. (eds) Biomedical Image Registration. WBIR 2024. Lecture Notes in Computer Science, vol 15249. Springer, Cham. [https://doi.org/10.1007/978-3-031-73480-9_25](https://doi.org/10.1007/978-3-031-73480-9_25)
