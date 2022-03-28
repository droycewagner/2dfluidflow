/*******************************************************************************

Full documentation is given in the implementation of these classes at mesh2d.cpp

struct rgb_tuple:
  A simple structure to represent a color as an rgb triple

struct Boundary:
  This is supposed to represent a boundary condition along one side of the mesh
  grid. The double 'value' represents the value of the boundary; the enum
  boundary_type allows the user to specify a Dirichlet or Neumann boundary
  condition.

class Differentiator:
  The Differentiator class stores differentials dx, dy as doubles and has member
  functions to compute any first or second partial derivatives of a matrix.

struct mesh2d:
  This is supposed to represent a 2d grid simulating a rectangular region
  containing fluid.

*******************************************************************************/

#ifndef MESH2D_H
#define MESH2D_H

#include <eigen3/Eigen/Core>

struct rgb_tuple {
public:
  rgb_tuple ();
  rgb_tuple (const float r, const float g, const float b);
  float r();
  float g();
  float b();
private:
  float red, green, blue;
};


struct Boundary {
public:
  enum boundary_type {Dirichlet, Neumann};
  Boundary(const boundary_type _type, const double _value);
  Boundary();
  double getValue();
  boundary_type getType();
private:
  double value;
  boundary_type type;
};


class Differentiator {
private:
  double dx, dy;
public:
  Differentiator();
  Differentiator(double _dx, double _dy);
  //use of templates allows one to take matrix derivatives of most Eigen types
  template <typename T> auto x(const Eigen::DenseBase<T>& m);
  template <typename T> auto y(const Eigen::DenseBase<T>& m);
  template <typename T> auto xx(const Eigen::DenseBase<T>& m);
  template <typename T> auto xy(const Eigen::DenseBase<T>& m);
  template <typename T> auto yy(const Eigen::DenseBase<T>& m);
};


struct mesh2d {
  using matrix=Eigen::MatrixXd;
private:
  matrix u, v, p;
  //having memory allocated for some derivatives seems to save time.
  matrix p_xy;
  matrix u_star, v_star;
  matrix u_star_x,v_star_y;
  double rho, mu, CFL;
  double tolerance;
  double dx, dy, dt;
  int nrow, ncol, nr, nc;
  Boundary u_left, u_right, u_top, u_bottom;
  Boundary v_left, v_right, v_top, v_bottom;
  Boundary p_left, p_right, p_top, p_bottom;
  Boundary default_boundary;
  Differentiator D;
  void resetUBoundary();
  void resetVBoundary();
  void resetPBoundary();
  void applyCFL();
public:
  mesh2d(const int _nrow, const int _ncol, const double length, const double breadth);
  void setFluid(const double _rho, const double _mu, const double _CFL);
  void setTolerance(const double _tolerance);
  double getDt();
  void setUBoundary(Boundary left, Boundary right, Boundary top, Boundary bottom);
  void setVBoundary(Boundary left, Boundary right, Boundary top, Boundary bottom);
  void setPBoundary(Boundary left, Boundary right, Boundary top, Boundary bottom);
  void writeFile(const std::string& filename);
  void writeImage(const std::string& filename, const double min, const double max);
  double maxVel();
  double minVel();
  void doIteration();
  matrix getU();
  matrix getV();
  matrix getP();
};

#endif
