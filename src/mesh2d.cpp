#include "mesh2d.h"
#include <chrono>
#include <Magick++.h>
#include <iostream>
#include <tuple>

//shorthand for a column/row of ones.
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

void mesh2d::set_CFL (const double CFL1) {
  CFL=CFL1;
}

void mesh2d::set_dt () {
  double denom=(u.array().maxCoeff()/dx+v.array().maxCoeff()/dy);
  dt=(denom==0)?CFL*(dx+dy):CFL/denom;
}

double mesh2d::get_dt() {return dt;}

//Saves vertical velocity boundary conditions to the mesh2d structure
void mesh2d::SetUBoundary(Boundary left, Boundary right, Boundary top, Boundary bottom) {
  u_l=left;u_r=right;u_t=top;u_b=bottom;
  ResetUBoundary(left,right,top,bottom);
}

//Set boundary conditions for horizontal velocity
void mesh2d::ResetUBoundary(Boundary left, Boundary right, Boundary top, Boundary bottom) {
  if (left.deriv()==0) u.col(0)=left.value()*onec(nrow);
    else u.col(0)=-left.value()*dx*onec(nrow)+u.col(1);
  if (right.deriv()==0) u.col(ncol-1)=right.value()*onec(nrow);
    else u.col(ncol-1)=right.value()*dx*onec(nrow)+u.col(ncol-2);
  if (top.deriv()==0) u.row(nrow-1)=2*top.value()*oner(ncol)-u.row(nrow-2);
    else u.row(nrow-1)=-top.value()*dy*oner(ncol)+u.row(nrow-2);
  if (bottom.deriv()==0) u.row(0)=2*bottom.value()*oner(ncol)-u.row(1);
    else u.row(0)=bottom.value()*dy*oner(ncol)+u.row(1);
}

//Saves vertical velocity boundary conditions to the mesh2d structure
void mesh2d::SetVBoundary(Boundary left, Boundary right, Boundary top, Boundary bottom) {
  v_l=left;v_r=right;v_t=top;v_b=bottom;
  ResetVBoundary(left,right,top,bottom);
}

//Set boundary conditions for vertical velocity
void mesh2d::ResetVBoundary(Boundary left, Boundary right, Boundary top, Boundary bottom) {
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

//Compute 1st derivative in y
matr mesh2d::D_y(const matr& m) {
  return (m.block(2,1,nr,nc)
    -m.block(0,1,nr,nc))/(2*dy);
}

//Compute 1st derivative in x
matr mesh2d::D_x(const matr& m) {
  return (m.block(1,2,nr,nc)
    -m.block(1,0,nr,nc))/(2*dx);
}

//Compute 2nd derivative in y
matr mesh2d::D_yy(const matr& m) {
  return (m.block(2,1,nr,nc)
    -2*m.block(1,1,nr,nc)
    +m.block(0,1,nr,nc))/(pow(dy,2));
}

//Compute 2nd derivative in x
matr mesh2d::D_xx(const matr& m) {
  return (m.block(1,2,nr,nc)
    -2*m.block(1,1,nr,nc)
    +m.block(1,0,nr,nc))/(pow(dx,2));
}

//Compute mixed partial derivative
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

//updates u_star and v_star in the mesh2d structure
void mesh2d::getStarredVelocities() {
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

//Solves the Poisson equation to get pressure matrix p
void mesh2d::SolvePoisson() {
  double error=1;
  double tol=.001;
  double factor=1/(2/pow(dx,2)+2/pow(dy,2));

  for (int i=0;i<=500;++i) {
    if (tol>=error) {
      break;}
    matr p_old=p;
    p_xy=D_xy(p);
    p.block(1,1,nr,nc)=p_xy*factor-
      (rho*factor/dt)*(D_x(u_star)+D_y(v_star));
    error=(p-p_old).array().abs().maxCoeff();
    ResetPBoundary(p_l,p_r,p_t,p_b);//necessary if there are Neumann conditions.
  }
}

//Last step in solving Navier-Stokes; solve the momentum equation;
//updates u,v in the mesh2d structure.
void mesh2d::SolveMomentum() {
  u.block(1,1,nr,nc)=u_star.block(1,1,nr,nc)-(dt/rho)*D_x(p);
  v.block(1,1,nr,nc)=v_star.block(1,1,nr,nc)-(dt/rho)*D_y(p);
}

void mesh2d::write2file(std::string filename) {
  std::ofstream outfile;
  outfile.open (filename);
  for (int r=1;r<nrow-1;++r) { for (int c=1;c<ncol-1;++c) {
    outfile<<p(r,c)<<"\t"<<u(r,c)<<"\t"<<v(r,c)<<"\n";
  }}
  outfile.close();
  return;
}

matr mesh2d::get_u() {return u;}

std::tuple<double,double,double> mesh2d::rainbow_scale(double vl, double min, double max) {
  double val=(vl-min)/(max-min)*5;//scales from 0-6 if min<vl<max
  val=5-val;
  int maxrgb=1;
  if (val<1) return {maxrgb,val*maxrgb,0};//red -> yellow
  else if (val<2) return {maxrgb*(2-val),maxrgb,0};//yellow -> green
  else if (val<3) return {0,maxrgb,maxrgb*(val-2)};//green -> cyan
  else if (val<4) return {0,maxrgb*(4-val),maxrgb};//cyan -> blue
  else if (val<=5) return {maxrgb*(val-4),0,maxrgb};//blue -> magenta
  //else if (val<=6) return {maxrgb,0,maxrgb*(6-val)};//magenta -> red
  else return {maxrgb,maxrgb,maxrgb};//white if val is outside usual range
}

void mesh2d::write2image(std::string filename, double min, double max) {

  //create a grid to place flow vectors
  std::vector<int> xpos, ypos;
  int len=7;int sp=2;
  for (int i=len+sp;i<nr-sp;i+=2*len+sp) xpos.push_back(i);
  for (int i=len+sp;i<nc-sp;i+=2*len+sp) ypos.push_back(i);

  //find speed at eaceh point of the grid
  matr vel=(u.cwiseProduct(u)
    +v.cwiseProduct(v))
    .array().sqrt().matrix().block(1,1,nr,nc);

  //create image
  Magick::Image image( Magick::Geometry(nr,nc),Magick::Color(MaxRGB,MaxRGB,MaxRGB,0));

  //create rainbow-colored background to represent velocity at each grid point.
  std::tuple<double,double,double> col;
  for (int i=0;i<nr;i++) {for (int j=0;j<nc;j++) {
    col=rainbow_scale(vel.coeff(i,j),min,max);
    image.pixelColor(i,j,Magick::Color(std::get<0>(col)*MaxRGB,std::get<1>(col)*MaxRGB,
      std::get<2>(col)*MaxRGB,MaxRGB));
  }}

  //add lines representing velocity vectors along a grid
  image.strokeColor("black");
  image.fillColor("black");
  image.strokeWidth(.5);
  std::list<Magick::Drawable> drawList;
  std::pair<double,double> vv;
  double vel_t;
  for (auto & x : xpos) {for (auto & y : ypos) {
    vv={u.coeff(x+1,y+1),v.coeff(x+1,y+1)};
    if (vv.first!=0&&vv.second!=0)
    drawList.push_back(Magick::DrawableLine(x,y,x+int(len*vv.first/vel.coeff(x,y)),y+int(len*vv.second/vel.coeff(x,y))));
  }}
  image.draw(drawList);

  //add a white dot at each point of the grid where we have placed a vector above
  for (auto & x : xpos) {for (auto & y : ypos) {
    image.pixelColor(x,y,Magick::Color(MaxRGB,MaxRGB,MaxRGB,MaxRGB));
  }}

  //write image to file
  try {
    image.write( filename );
  } catch( std::exception &error_ ) {
    std::cout << "Couldn't write to file "+filename << error_.what() << std::endl;
    return;
  }
}

double mesh2d::max_vel() {
return (u.block(1,1,nr,nc).cwiseProduct(u.block(1,1,nr,nc))
  +v.block(1,1,nr,nc).cwiseProduct(v.block(1,1,nr,nc)))
  .array().sqrt().maxCoeff();
}

double mesh2d::min_vel() {
return (u.block(1,1,nr,nc).cwiseProduct(u.block(1,1,nr,nc))
  +v.block(1,1,nr,nc).cwiseProduct(v.block(1,1,nr,nc)))
  .array().sqrt().minCoeff();
}

//performs an iteration of one time step
void mesh2d::do_iteration() {
  set_dt();
  //setting the boundary during each iteration is necessary if there are Neumann conditions.
  ResetUBoundary(u_l,u_r,u_t,u_b);
  ResetVBoundary(v_l,v_r,v_t,v_b);
  ResetPBoundary(p_l,p_r,p_t,p_b);
  getStarredVelocities();
  SolvePoisson();
  SolveMomentum();
}
