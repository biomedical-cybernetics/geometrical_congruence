### Description

The folder "patent_geometry" contains the main function "compute_GC_GRE.m" to computed geometrical congruence (GC) and greedy routing efficiency (GRE) in networks with patent geometry. The example script "run_example.m" loads nPSO example networks stored in "nPSO_example_data.mat" and plots results for a case scenario.
The folder "latent_geometry_geometry" contains the main function "compute_GC_latent.m" to computed geometrical congruence (GC) in networks with latent geometry. The example script "run_example.m" loads the brain connectomes dataset stored in "brain_connectomes_data.mat" and plots results for gender and age analyses.
The folder "plot_hyperbolic_network" contains the main function "plot_hyperbolic_network.m" to plot networks in the hyperbolic space (native or poincare unit disk representations). The example script "run_example_figure_1.m" loads a nPSO example network stored in "nPSO_example_data.mat" and reproduces two panels of figure 1.
See the documentation at the beginning of the functions for input and output details.
See the documentation at the beginning of the example scripts for more details on the data and analyses.
The source files with ".c" extension needs to be compiled to build MEX support functions, more details in the section below.


### Build MEX functions

MEX functions for Windows (".mexw64" extension) and Linux (".mexa64" extension) are already provided.
If they are not working, please build them again from the C source files.

- Windows
Go to MATLAB "Add-Ons" and install "MATLAB Support for MinGW-w64 C/C++ Compiler"
Build the MEX functions using the following MATLAB commands (change the MinGW path if needed):
mex C:\ProgramData\MATLAB\SupportPackages\R2020b\3P.instrset\mingw_w64.instrset\lib\gcc\x86_64-w64-mingw32\6.3.0\libgomp.a compute_pGRP_mex.c CFLAGS='$CFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp'
mex C:\ProgramData\MATLAB\SupportPackages\R2020b\3P.instrset\mingw_w64.instrset\lib\gcc\x86_64-w64-mingw32\6.3.0\libgomp.a compute_pTSP_mex.c CFLAGS='$CFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp'
mex meanratio2_mask_mex.c
It will generate the following files: compute_pGRP_mex.mexw64, compute_pTSP_mex.mexw64, meanratio2_mask_mex.mexw64.

- Linux
Build the MEX functions using the following MATLAB commands:
mex compute_pGRP_mex.c CFLAGS='$CFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp'
mex compute_pTSP_mex.c CFLAGS='$CFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp'
mex meanratio2_mask_mex.c
It will generate the following files: compute_pGRP_mex.mexa64, compute_pTSP_mex.mexa64, meanratio2_mask_mex.mexa64.

- Mac
Building the MEX functions has not been tested yet, a procedure similar to Linux should work.
It will generate the following files: compute_pGRP_mex.mexmaci64, compute_pTSP_mex.mexmaci64, meanratio2_mask_mex.mexmaci64.


### Reference

"Geometrical congruence, greedy navigability and myopic transfer in complex networks and brain connectomes"
C. V. Cannistraci, A. Muscoloni, Nature Communications, 2022 (accepted for publication)


### Contact

Alessandro Muscoloni: alessandro.muscoloni@gmail.com
Carlo Vittorio Cannistraci: kalokagathos.agon@gmail.com
