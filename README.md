# 2dfluidflow

Two-Dimensional Laminar Fluid Flow
==================================
downloaded from [DRW Github](github.com)

This program gives a numerical solution and visualization (TODO) for the 2d Navier-Stokes equations in single phase laminar flow. A "mesh2d" data type is provided (see _src/mesh2d.cpp_, _src/mesh2d.h_) which allows users to set a variety of initial conditions on a rectangle of chosen shape and resolution, as well as set the properties of the fluid. A sample implementation is given in _fluidflow.cpp_ where a lid cavity test is performed, meeting a benchmark of Ghia et al. (1982) to a relative error of about 0.1%.

If the system requirements are met (see the appropriate section below) this can be compiled/run by navigating to the root directory of the project and running

    make
    ./fluidflow

This is a C++ adaptation of Gaurav Deshmukh's project available on
[Github](https://github.com/gauravsdeshmukh/FlowPy)
with some exposition at
[towardsdatascience.com](https://towardsdatascience.com/computational-fluid-dynamics-using-python-modeling-laminar-flow-272dad1ebec).


Differences from Predecessor
----------------------------
1. When no graphics are output, this implementation runs in just under 1/3 the time of the above reference using the 150 second 257x257 lid cavity test. This puts operation speeds well under real-time.
2. Use of C++ classes and the elimination of code duplication in the computation of derivatives reduces program length by more than one third.
* Error fixed in SetPBoundary (see _mesh2d.h_), which happens not to effect Gaurav's lid cavity test due to its particular initial conditions.
* A custom visualization is used in place of python's mathplotlib (TODO).


System Requirements
-------------------
I have both compiled and run this on my own machine running Ubuntu 20.04.

This codes uses [Eigen](https://eigen.tuxfamily.org/) to implement matrices. Installation procedures and documentation can be found [here](https://eigen.tuxfamily.org/dox/GettingStarted.html).

The use of std::filesystem requires C++17 or later during compilation.

__Important__: the use of -Ofast (or -O3) flag in the makefile is crucial; this gives a runtime speedup on the order of 20x in the use of Eigen operations.

If you are compiling this, the final visualization requires imagemagick/magick++. On some linux systems, you can meet this dependency by installing the packages:

    sudo apt install imagemagick libgraphicsmagick++1-dev graphicsmagick-libmagick-dev-compat


TODO
----
1. This code makes no guarantee of exception safety.
2. There is an unfulfilled promise of graphics output -- coming soon!
