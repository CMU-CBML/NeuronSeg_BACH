NEUROSEG_BACH Software Tool: This software performs 2D and 3D automatic segmentation of neuron cell fluorescence microscopy images.  We have divided 2D and 3D image segmentation into separate folders. The modules for the software in each folder are given as follows:

BsplineComposition: Functions related to composition update and smooth representation of level set function using B-splines is provided.

Initialization: In this folder the code for automatic initilization is provided. For the images given in the dates, corresponding initialization functions are provided.

IterationLoopFunctions: The functions used in each iteration update re provided in this folder.

LevelSetFunctions: The functions used in the local image fitting energy functional are provided folder.

SetParameters: In this folder, the user-defined parameters are set in the code associated with each dataset.

THBSpline: In this folder, functions associated with THBspline refinement are provided in this folder.

data: This folder has all the images used in the demo code.

The code is tested with the MATLAB 2019a version. For the slower functions, the code is compiled using C++ functions and run using MEX functions. 

General instructions to run the code:
1. Input image that is to be segmented is to be save in data/ folder.
2. To run the demo code, run the files starting with "mainfunction" for each numerical example shown in paper.
3. Before running the mainfunction.m file for 3D images, first run compilemex.m program to compile the MEX functions to increase the speed. The, run the mainfunction script. This script compiles all the MEX functions using the OpenMP parallel framework, provided that a supported compiler is installed. If compilation fails, the native Matlab .m files can be used by removing the "_mex" suffix from the function names.

For different examples, we have used different main files for running with the respective settings:

1. 2D neuron fluorescence microscopy images examples (Fig. 2 in article)
run file: mainfunction_neuron01.m (01,02,03,04) script is used to run these examples. 

The images used are obtained from the software code: https://github.com/kilho/NIA (01.tif, 06.tif, 08.tif, 12.tif)

The image files should be stored in data/ folder. The intensities are adjusted to be in the range [0, 255]. User can use the function neuron01_initilization.m function. The user-defined parameters are set in the file setparameters_neuron01_seg.m script in the SetParameters folder. 


2.Neuron Synapse fluorescence microscopy examples (Fig. 3 & 4 in article):
run file: mainfunction_neuronsynapse_1.m, mainfunction_neuronsynapse_2.m scripts are used to run these examples. 

The images used are obtained from the link: https://doi.org/doi:10.7295/W9CIL36160 and https://doi.org/doi:10.7295/W9CIL12435

The image files should be stored in data/ folder. The intensities are adjusted to be in the range [0, 255]. User can use the function neuronsynapse_1_initilization.m function. The user-defined parameters are set in the file setparameters_neuronsynapse_1_seg.m script in the SetParameters folder.  

3. 3D Synthetic example of Helix geometry (Fig. 5 in article)
run file: mainfunction_helix.m script is used to run this example.

The image file 'helix_img.mat' is stored in data/ folder. 

The intensities are adjusted to be in the range [0, 255]. The user-defined parameters are set in the file setparameters_helix.m script in the SetParameters folder.  


4. 3D Neuron Stacked Image example of olfactory projection fiber image data obtained from the DIADEM challenge dataset (Figs. 7-8 in article)

run file: mainfunction_neuron_olfactoryfiber.m script is used to run this example.

The images used are obtained from the link: http://diademchallenge.org/olfactory_projection_fibers_readme.html

The image files should be stored in data/ folder. The intensities are adjusted to be in the range [0, 255]. User can use the function neuron_initilization.m function. The user-defined parameters are set in the file setparameters_olfactoryfiber.m script in the SetParameters folder. 

To run this example, one needs to download the TREES Toolbox and save it in the folder OlfactoryFiber/ folder. The TREES toolbox is to be downloaded from: https://www.treestoolbox.org/download.html

The gold standard reconstructions are in .SWC format and can be downloaded from: http://diademchallenge.org/olfactory_projection_fibers_readme.html
Extract 'Olfactory Projection Fibers.rar'. The gold standard standard reconstructions can be found in the subfolder: 'Gold Standard Reconstructions'. Save the files in OlfactoryFiber/ folder.



