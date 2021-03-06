# 2dfluidflow

Two-Dimensional Laminar Fluid Flow
==================================

This program gives a numerical solution and visualization for the 2d Navier-Stokes equations in single phase laminar flow. A "Mesh2D" data type is provided (see _src/mesh2d.cpp_, _src/mesh2d.h_) which allows users to set a variety of initial conditions on a rectangle of chosen shape and resolution, as well as set the properties of the fluid. A sample simulation is given in _fluidflow.cpp_ where a lid cavity test is performed, meeting a benchmark of Ghia et al. (1982) to a relative error of about 6%.

If the system requirements are met (see the appropriate section below) this can be compiled/run by navigating to the root directory of the project and running

    make
    ./fluidflow

Detailed documentation is given in the implementation _src/mesh2d.cpp_.


Output Options
--------------

A method _writeFile()_ is provided to print tab-delimited velocity and pressure matrices to file and there is a method _writeImage_ to produce colored images representing the fluid flow. The following options can be set for visual output in the file _fluidflow.cpp_:

* _bool write_file_: when set to 1, writes tab-delimited pressure/velocity matrices to file in a folder named _res_.
* _bool write_image_: when set to 1, produces png files in the _res_ folder, representing the fluid flow
* _double fps_: this will indicate how many times per second (in simulation time) an image/file should be produced
* _double min_col, max_col_: these give the coloration scale. Images contain colors from magenta (representing zero velocity) to red (representing max velocity). A velocity below this range will be colored white and a velocity above the maximum will be colored black.

The function _writeImage_ does the following:

1. creates an image with the same dimensions as the Mesh2D object,
2. colors each pixel by speed, based on _rainbow_scale_,
3. places short black lines (vectors) to represent the velocity direction (if the velocity is nonzero),
4. places a white dot at the base of each vector.

For example, the following image is produced from the lid cavity test near the 25 second mark.

![img_050](https://user-images.githubusercontent.com/65519600/163690494-c56e42ad-43ef-44fd-b89a-d16f9762a505.png)


Once images have been produced, ffmpeg users might create a lossless video from the output:

    ffmpeg -framerate 20 -pattern_type glob -i "res/*.png" -c:v libx264 -crf 0 output.mp4


Differences from Predecessor
----------------------------
This is originally a C++ adaptation of Gaurav Deshmukh's project available on
[Github](https://github.com/gauravsdeshmukh/FlowPy)
with some exposition at
[towardsdatascience.com](https://towardsdatascience.com/computational-fluid-dynamics-using-python-modeling-laminar-flow-272dad1ebec).

* When no graphics are output, this implementation runs in under 1/6-th the time of the above reference using the 150 second 257x257 lid cavity test. This puts operation speeds well under real-time. The graphics generation also runs under real-time.
* Concision is gained from the use of C++ classes and the elimination of code duplication.
* Error fixed in SetPBoundary (see _mesh2d.h_), which happens not to effect the lid cavity test due to its particular initial conditions.
* A custom visualization is used in place of python's mathplotlib.


System Requirements
-------------------
I have both compiled and run this on my own machine running Ubuntu 20.04.

This code uses [Eigen](https://eigen.tuxfamily.org/) to implement matrices. Installation procedures and documentation can be found [here](https://eigen.tuxfamily.org/dox/GettingStarted.html).

The use of std::filesystem requires C++17 or later during compilation.

__Important__: use of the -Ofast (or -O3) flag in the makefile is crucial; this gives a runtime speedup on the order of 20x in the use of Eigen operations.

The final visualization requires imagemagick/magick++. On some Linux systems, you can meet this dependency by installing the packages:

    sudo apt install imagemagick libgraphicsmagick++1-dev graphicsmagick-libmagick-dev-compat
