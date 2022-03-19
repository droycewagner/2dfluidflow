/**********
Benchmark Lid Cavity Test:

In this simulation, a square 2d-fluid measuring 4x4 in length/breadth is
submitted to a flow of speed 1 along one side. We run the simulation, produce
graphical output and check the numerical result against the benchmark given by
Ghia et al. in

  Ghia, U.K.N.G., Ghia, K.N. and Shin, C.T., 1982.
  High-Re solutions for incompressible flow using the Navier-Stokes equations
  and a multigrid method. Journal of computational physics, 48(3), pp.387-411.

This should run in ~1 min on a reasonably fast machine (no output).
Final output for the benchmark should give
MSE=0.000268503
Average relative error: .00098211
**********/

/*******************
to convert png files to a lossless mp4:
ffmpeg -framerate 20 -pattern_type glob -i "res/*.png" -c:v libx264 -crf 0 output.mp4
*******************/

#include <cmath>
#include <iostream>
#include <filesystem>
#include <chrono>
#include <algorithm>
#include <Magick++.h>
#include "src/mesh2d.h"

bool make_file=0;//write pressure and velocity matrices to file?
bool make_vid=0;//make a video?
double fps=2;//approx frame rate for video output
double min_col=0;//sets color scale for image output
double max_col=0.99;

//set constants for dimension, runtime, fluid properties
double length=4;
double breadth=4;
double colpts=257;
double rowpts=257;
double sim_time=150;
double CFL=0.8;
double rho=1;
double mu=0.01;

//create boundary conditions
int Dirichlet=0;
int Neumann=1;
Boundary flow=Boundary(Dirichlet,1);
Boundary noslip=Boundary(Dirichlet,0);
Boundary zeroflux=Boundary(Neumann,0);
Boundary pressureatm=Boundary(Dirichlet,0);

// given int n,m returns a string which is the base-10 expression of
// n with zeros prepended, to be of length m.
std::string pad_int(int n, int m) {
  std::string str1=std::to_string(n);
  std::string str2=(str1.length()>=m)?"":std::string(m-str1.length(),'0');
  return str2+str1;
}


int main() {
  //gives Magick package knowledge of the working directory.
  Magick::InitializeMagick(nullptr);
  //std::ios::sync_with_stdio(false);

  //initialize the mesh2d object
  mesh2d mesh=mesh2d(rowpts,colpts);
  mesh.setDims(4,4);
  mesh.setFluid(rho,mu);
  mesh.set_CFL(CFL);
  mesh.SetUBoundary(noslip,noslip,flow,noslip);
  mesh.SetVBoundary(noslip,noslip,noslip,noslip);
  mesh.SetPBoundary(zeroflux,zeroflux,pressureatm,zeroflux);

  //create an empty directory /res for output
  if (make_vid||make_file) {
    try {
      if (!std::filesystem::is_directory("res") || !std::filesystem::exists("res"))
        std::filesystem::create_directory("res");
      for (auto& path: std::filesystem::directory_iterator("res"))
        std::filesystem::remove_all(path);
    } catch (const std::exception& e) {
       std::cout << e.what();
       std::cout<<"couldn't make output directory -- won't produce output\n";
       make_vid=0;make_file=0;
    }
  }

  auto start = std::chrono::high_resolution_clock::now();
  double t=0;int count=0;int padn=std::to_string(int(sim_time*fps+1)).length();

  //run the simulation
  while (t<sim_time) {
    mesh.do_iteration();
    if (t>=(1/fps)*count) {
      count++;
      std::cout<<"\rSimulation time is "<<t<<"       "<<std::flush;
      if (make_vid||make_file) {
        if (make_vid)
          mesh.write2image("res/img_"+pad_int(count,padn)+".png",min_col,max_col);
        if (make_file)
          mesh.write2file("res/puv_"+std::to_string(count));
      }
    }
    t+=mesh.get_dt();
  }

  std::cout<<"\n";
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration<float>(stop - start);
  std::cout <<"Time elapsed: "<< duration.count()<<"sec\n\n";

  /************************
  perform benchmark test:
  ************************/
  //Ghia's 16 velocities with corresponding positions along the vertical median
  std::vector<double> x_g,y_g,u_g,v_g;
  y_g={0,0.0547,0.0625,0.0703,0.1016,0.1719,0.2813,0.4531,0.5,0.6172,0.7344,0.8516,0.9531,0.9609,0.9688,0.9766};
  u_g={0,-0.08186,-0.09266,-0.10338,-0.14612,-0.24299,-0.32726,-0.17119,-0.11477,0.02135,0.16256,0.29093,0.55892,0.61756,0.68439,0.75837};
  x_g={0,0.0625,0.0703,0.0781,0.0983,0.1563,0.2266,0.2344,0.5,0.8047,0.8594,0.9063,0.9453,0.9531,0.9609,0.9688};
  v_g={0,0.1836,0.19713,0.20920,0.22965,0.28124,0.30203,0.30174,0.05186,-0.38598,-0.44993,-0.23827,-0.22847,-0.19254,-0.15663,-0.12146};

  //compute error from benchmark
  double mse=0;
  double rel=0;
  auto mid_u=mesh.get_u().block(1,floor(colpts/2),rowpts,1).array();
  for (int i=0;i<u_g.size();i++) {
    mse+=pow(u_g[i]-mid_u.coeff(floor(y_g[i]*rowpts),0),2);
    rel+=abs(u_g[i]);
  }

  //print error from benchmark
  std::cout <<"Ghia et al. Cavity test benchmark results:\n";
  std::cout<<"  MSE for horizontal component of velocity along vertical median: "
    <<mse/u_g.size()<<"\n";
  std::cout<<"  Average relative error: "
    <<mse/rel<<"\n\n";

  return 0;
}
