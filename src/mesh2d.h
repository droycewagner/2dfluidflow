/*******************************************************************************

Full documentation is given in the implementation of these classes at mesh2d.cpp

struct RGBTuple:
  A simple structure to represent a color as an rgb triple

struct Boundary:
  This is supposed to represent a boundary condition along one side of the mesh
  grid. The double 'value' represents the value of the boundary; the enum class
  Condition allows the user to specify a Dirichlet or Neumann boundary
  condition.

struct Mesh2D:
  This is supposed to represent a 2d grid simulating a rectangular region
  containing fluid.

*******************************************************************************/

#ifndef MESH2D_H
#define MESH2D_H

#include <eigen3/Eigen/Core>

struct RGBTuple {
public:
  RGBTuple ();
  RGBTuple (const float r, const float g, const float b);
  float r() const;
  float g() const;
  float b() const;
private:
  float red=0;
  float green=0;
  float blue=0;
};


struct Boundary {
public:
  enum class Condition {Dirichlet, Neumann};
  Boundary(const Condition _type, const double _value);
  Boundary();
  double getValue() const;
  Condition getType() const;
private:
  double value;
  Condition type;
};


struct Mesh2D {
  using matrix=Eigen::MatrixXd;
private:
  matrix u, v, p;
  matrix p_xy;
  matrix u_star, v_star;
  matrix u_star_x,v_star_y;
  double rho, mu, CFL;
  double tolerance;
  double dx, dy, dt;
  int nrow, ncol;
  Boundary u_left, u_right, u_top, u_bottom;
  Boundary v_left, v_right, v_top, v_bottom;
  Boundary p_left, p_right, p_top, p_bottom;
  void resetUBoundary();
  void resetVBoundary();
  void resetPBoundary();
public:
  Mesh2D(const int _nrow, const int _ncol, const double length, const double breadth);
  void setFluid(const double _rho, const double _mu, const double _CFL);
  void setTolerance(const double _tolerance);
  void setUBoundary(const Boundary left, const Boundary right, const Boundary top, const Boundary bottom);
  void setVBoundary(const Boundary left, const Boundary right, const Boundary top, const Boundary bottom);
  void setPBoundary(const Boundary left, const Boundary right, const Boundary top, const Boundary bottom);
  void writeFile(const std::string& filename) const;
  void writeImage(const std::string& filename, const double min, const double max) const;
  void doIteration();
  double getDt() const;
  const matrix& getU() const;
  const matrix& getV() const;
  const matrix& getP() const;
};

#endif
