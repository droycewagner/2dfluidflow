/*******************************************************************************
Description of Boundary class:
  This is supposed to represent a boundary condition along one side of the mesh
  grid. The double val represents the value of the boundary; The int deriv is
  supposed to indicate the derivate to which the condition applies; we use
  deriv=0 to represent a Dirichlet boundary condition and deriv=1 to represent
  a Neumann boundary condition.

Description of mesh2d struct:
At its core, this is supposed to represent a 2d grid, which in turn represents
a rectangular region containing fluid.
  u: a matrix giving horizontal velocity at each point of the grid
  v: same, for vertical velocity
  p: same, for pressure
  rho, mu, CFL: these are constants describing the fluid
  dx,dy: these give the space between points of the mesh
  dt: the timestep
The various functions defined here are described further in mesh2d.cpp where 
they are implemented.

*******************************************************************************/

#ifndef MESH2D_H
#define MESH2D_H

#include <vector>
#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Core>

typedef Eigen::MatrixXd matr;

struct Boundary {
private:
  double val;
  int der;// val applies to "deriv"th derivative
public:
  Boundary (int der1,double val1) ;
  Boundary ();
  double value();
  int deriv();
};

struct mesh2d {
private:
  matr u, v, p, u_face, v_face, u_star, v_star;
  double rho, mu, CFL;
  double dx, dy, dt;
  int nrow, ncol, nr, nc;
  Boundary u_l, u_r, u_t, u_b;
  Boundary v_l, v_r, v_t, v_b;
  Boundary p_l, p_r, p_t, p_b;
  matr D_y(const matr& m);
  matr D_x(const matr& m);
  matr D_yy(const matr& m);
  matr D_xx(const matr& m);
  matr D_xy(const matr& m);
  void ResetUBoundary(Boundary left, Boundary right, Boundary top, Boundary bottom);
  void ResetVBoundary(Boundary left, Boundary right, Boundary top, Boundary bottom);
  void ResetPBoundary(Boundary left, Boundary right, Boundary top, Boundary bottom);
  void copyBoundary(const matr&a, matr&b);
  std::tuple<double,double,double> rainbow_scale(double vl, double min, double max);
  void getStarredVelocities();
  void SolvePoisson();
  void SolveMomentum();
  void set_dt();

public:
  mesh2d(int a, int b);
  void setDims(double length, double breadth);
  void setFluid(double rho1, double mu1);
  void set_CFL(const double CFL1);
  double get_dt();
  void SetUBoundary(Boundary left, Boundary right, Boundary bottom, Boundary top);
  void SetVBoundary(Boundary left, Boundary right, Boundary bottom, Boundary top);
  void SetPBoundary(Boundary left, Boundary right, Boundary bottom, Boundary top);
  matr get_u();
  void write2file(std::string filename);
  void write2image(std::string filename, double min, double max);
  double max_vel();
  double min_vel();
  void do_iteration();
};

#endif
