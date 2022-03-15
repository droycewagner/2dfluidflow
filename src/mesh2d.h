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
  matr u, v, p, u_face, v_face, u_star, v_star,
    u_x, u_y, u_xx, u_yy, v_x, v_y, v_xx, v_yy, u_star_x, v_star_y,
    p_x, p_y, p_xy;
  double rho, mu, CFL;
  double dx, dy, dt;
  int nrow, ncol, nr, nc;
  Boundary p_l, p_r, p_t, p_b;
  matr D_y(const matr& m);
  matr D_x(const matr& m);
  matr D_yy(const matr& m);
  matr D_xx(const matr& m);
  matr D_xy(const matr& m);
  void ResetPBoundary(Boundary left, Boundary right, Boundary top, Boundary bottom);
  void copyBoundary(const matr&a, matr&b);
public:
  mesh2d(int a, int b);
  void setDims(double length, double breadth);
  void setFluid(double rho1, double mu1);
  void set_dt(const double CFL);
  double get_dt();
  void SetUBoundary(Boundary left, Boundary right, Boundary bottom, Boundary top);
  void SetVBoundary(Boundary left, Boundary right, Boundary bottom, Boundary top);
  void SetPBoundary(Boundary left, Boundary right, Boundary bottom, Boundary top);
  void getStarredVelocities();
  void SolvePoisson();
  void SolveMomentum();
  matr get_u();
  void write2file(std::string filename);
};

#endif
