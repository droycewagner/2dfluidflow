#include "mesh2d.h"
#include <chrono>

//std::cout<<"varname\n"<<varname<<"\n"<<std::flush;

matr onec(int n) {return matr::Ones(n,1);}
matr oner(int n) {return matr::Ones(1,n);}

Boundary::Boundary(int der1,double val1) {
  val=val1;
  der=der1;
}

Boundary::Boundary () {}

double Boundary::value() {return val;}

int Boundary::deriv() {return der;}

mesh2d::mesh2d(int a, int b) {
  double rho, mu, CFL;
  nrow=a+2;ncol=b+2;
  nr=a;nc=b;
  u=matr::Zero(nrow,ncol);
  v=matr::Zero(nrow,ncol);
  p=matr::Zero(nrow,ncol);
  u_star=matr::Zero(nrow,ncol);
  v_star=matr::Zero(nrow,ncol);
  u_face=matr::Zero(a,b);
  v_face=matr::Zero(a,b);
}

void mesh2d::setDims (double length, double breadth) {
  dx=length/(nc-1);
  dy=breadth/(nr-1);
}

void mesh2d::setFluid(double rho1, double mu1) {
  rho=rho1;mu=mu1;
}

void mesh2d::set_dt (const double CFL) {
  double denom=(u.array().maxCoeff()/dx+v.array().maxCoeff()/dy);
  dt=(denom==0)?CFL*(dx+dy):CFL/denom;
}

double mesh2d::get_dt() {return dt;}

//Set boundary conditions for horizontal velocity
void mesh2d::SetUBoundary(Boundary left, Boundary right, Boundary top, Boundary bottom) {
  if (left.deriv()==0) u.col(0)=left.value()*onec(nrow);
    else u.col(0)=-left.value()*dx*onec(nrow)+u.col(1);
  if (right.deriv()==0) u.col(ncol-1)=right.value()*onec(nrow);
    else u.col(ncol-1)=right.value()*dx*onec(nrow)+u.col(ncol-2);
  if (top.deriv()==0) u.row(nrow-1)=2*top.value()*oner(ncol)-u.row(nrow-2);
    else u.row(nrow-1)=-top.value()*dy*oner(ncol)+u.row(nrow-2);
  if (bottom.deriv()==0) u.row(0)=2*bottom.value()*oner(ncol)-u.row(1);
    else u.row(0)=bottom.value()*dy*oner(ncol)+u.row(1);
}

//Set boundary conditions for vertical velocity
void mesh2d::SetVBoundary(Boundary left, Boundary right, Boundary top, Boundary bottom) {
  if (left.deriv()==0) v.col(0)=2*left.value()*onec(nrow)-v.col(1);
    else v.col(0)=-left.value()*dx*onec(nrow)+v.col(1);
  if (right.deriv()==0) v.col(ncol-1)=2*right.value()*onec(nrow)-v.col(ncol-2);
    else v.col(ncol-1)=right.value()*dx*onec(nrow)+v.col(ncol-2);
  if (top.deriv()==0) v.row(nrow-1)=top.value()*oner(ncol);
    else v.row(nrow-1)=-top.value()*dy*oner(ncol)+v.row(nrow-2);
  if (bottom.deriv()==0) v.row(0)=bottom.value()*oner(ncol);
    else v.row(0)=bottom.value()*dy*oner(ncol)+v.row(1);
}

//Saves pressure boundary conditions to the mesh2d structure for use in SolvePoisson(),
//then sets the boundary condition
void mesh2d::SetPBoundary(Boundary left, Boundary right, Boundary top, Boundary bottom) {
  p_l=left;p_r=right;p_t=top;p_b=bottom;
  ResetPBoundary(left,right,top,bottom);
}

//Set boundary conditions for pressure
void mesh2d::ResetPBoundary(Boundary left, Boundary right, Boundary top, Boundary bottom) {
  if (left.deriv()==0) p.col(0)=left.value()*onec(nrow);
    else p.col(0)=-left.value()*dx*onec(nrow)+p.col(1);
  if (right.deriv()==0) p.col(ncol-1)=right.value()*onec(nrow);//here's the error
    else p.col(ncol-1)=right.value()*dx*onec(nrow)+p.col(ncol-2);
  if (top.deriv()==0) p.row(nrow-1)=top.value()*oner(ncol);
    else p.row(nrow-1)=-top.value()*dy*oner(ncol)+p.row(nrow-2);
  if (bottom.deriv()==0) p.row(0)=bottom.value()*oner(ncol);
    else p.row(0)=bottom.value()*dy*oner(ncol)+p.row(1);
}

matr mesh2d::D_y(const matr& m) {
  return (m.block(2,1,nr,nc)
    -m.block(0,1,nr,nc))/(2*dy);
}

matr mesh2d::D_x(const matr& m) {
  return (m.block(1,2,nr,nc)
    -m.block(1,0,nr,nc))/(2*dx);
}

matr mesh2d::D_yy(const matr& m) {
  return (m.block(2,1,nr,nc)
    -2*m.block(1,1,nr,nc)
    +m.block(0,1,nr,nc))/(pow(dy,2));
}

matr mesh2d::D_xx(const matr& m) {
  return (m.block(1,2,nr,nc)
    -2*m.block(1,1,nr,nc)
    +m.block(1,0,nr,nc))/(pow(dx,2));
}

matr mesh2d::D_xy(const matr& m) {
  return (m.block(2,1,nr,nc)+m.block(0,1,nr,nc))/pow(dy,2)
    +(m.block(1,2,nr,nc)+m.block(1,0,nr,nc))/pow(dx,2);
}

//copies first/last row/column from a to b.
void mesh2d::copyBoundary(const matr&a, matr&b) {
  b.row(0)=a.row(0);
  b.col(0)=a.col(0);
  b.row(b.rows()-1)=a.row(a.rows()-1);
  b.col(b.cols()-1)=a.col(a.rows()-1);
  return;
}

void mesh2d::getStarredVelocities() {
  typedef std::chrono::high_resolution_clock Time;
  typedef std::chrono::duration<float> fsec;
  u_face=(u.block(1,1,nr,nc)+u.block(1,2,nr,nc)
    +u.block(0,1,nr,nc)+u.block(0,2,nr,nc))/4;
  v_face=(v.block(1,1,nr,nc)+v.block(1,0,nr,nc)
    +v.block(2,1,nr,nc)+v.block(2,0,nr,nc))/4;
  copyBoundary(u,u_star);
  copyBoundary(v,v_star);
  u_star.block(1,1,nr,nc)=u.block(1,1,nr,nc)
    -dt*(u.block(1,1,nr,nc).cwiseProduct(D_x(u))+v_face.cwiseProduct(D_y(u)))
    +dt*(mu/rho)*(D_xx(u)+D_yy(u));
  v_star.block(1,1,nr,nc)=v.block(1,1,nr,nc)
    -dt*(u_face.cwiseProduct(D_x(v))+v.block(1,1,nr,nc).cwiseProduct(D_y(v)))
    +dt*(mu/rho)*(D_xx(v)+D_yy(v));
}

void mesh2d::SolvePoisson() {
  double error=1;
  double tol=.001;
  double factor=1/(2/pow(dx,2)+2/pow(dy,2));

  for (int i=0;i<=500;++i) {
    if (tol>=error) {//std::cout<<"went "<<i<<" times.\n"<<std::flush;
      break;}
    matr p_old=p;
    p_xy=D_xy(p);
    p.block(1,1,nr,nc)=p_xy*factor-
      (rho*factor/dt)*(D_x(u_star)+D_y(v_star));
    error=(p-p_old).array().abs().maxCoeff();
    ResetPBoundary(p_l,p_r,p_t,p_b);//necessary if there are Neumann conditions.
    //copyBoundary(p_old,p);
  }
}

void mesh2d::SolveMomentum() {
  u.block(1,1,nr,nc)=u_star.block(1,1,nr,nc)-(dt/rho)*D_x(p);
  v.block(1,1,nr,nc)=v_star.block(1,1,nr,nc)-(dt/rho)*D_y(p);
}

matr mesh2d::get_u() {return u;}

void mesh2d::write2file(std::string filename) {
  std::ofstream outfile;
  outfile.open (filename);
  for (int r=1;r<nrow-1;++r) { for (int c=1;c<ncol-1;++c) {
    outfile<<p(r,c)<<"\t"<<u(r,c)<<"\t"<<v(r,c)<<"\n";
  }}
  outfile.close();
  return;
}
