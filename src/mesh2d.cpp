#include "mesh2d.h"
#include <Magick++.h>
#include <iostream>
#include <fstream>

rgb_tuple::rgb_tuple(): red(0),green(0),blue(0) {}

rgb_tuple::rgb_tuple(const float _r, const float _g, const float _b) {
  red=_r;
  green=_g;
  blue=_b;
}

float rgb_tuple::r() {return red;}

float rgb_tuple::g() {return green;}

float rgb_tuple::b() {return blue;}


Boundary::Boundary(boundary_type _type, double _value) {
  value=_value;
  type=_type;
}

Boundary::Boundary(): type(Boundary::Dirichlet), value(0) {}

Boundary::boundary_type Boundary::getType() {return type;}

double Boundary::getValue() {return value;}


Differentiator::Differentiator(): dx(1), dy(1) {}

Differentiator::Differentiator(double _dx, double _dy) {
  dx=_dx;
  dy=_dy;
}

//Compute 1st derivative in y
template <typename T> auto Differentiator::y(const Eigen::DenseBase<T>& m) {
  return (m.block(2,1,m.rows()-2,m.cols()-2)-m.block(0,1,m.rows()-2,m.cols()-2))/(2*dy);
}

//Compute 1st derivative in x
template <typename T> auto Differentiator::x(const Eigen::DenseBase<T>& m) {
  return (m.block(1,2,m.rows()-2,m.cols()-2)-m.block(1,0,m.rows()-2,m.cols()-2))/(2*dx);
}

//Compute 2nd derivative in y
template <typename T> auto Differentiator::yy(const Eigen::DenseBase<T>& m) {
  return (m.block(2,1,m.rows()-2,m.cols()-2)-2*m.block(1,1,m.rows()-2,m.cols()-2)
    +m.block(0,1,m.rows()-2,m.cols()-2))/(pow(dy,2));
}

//Compute 2nd derivative in x
template <typename T> auto Differentiator::xx(const Eigen::DenseBase<T>& m) {
  return (m.block(1,2,m.rows()-2,m.cols()-2)-2*m.block(1,1,m.rows()-2,m.cols()-2)
    +m.block(1,0,m.rows()-2,m.cols()-2))/(pow(dx,2));
}

//Compute mixed partial derivative
template <typename T> auto Differentiator::xy(const Eigen::DenseBase<T>& m) {
  return (m.block(2,1,m.rows()-2,m.cols()-2)+m.block(0,1,m.rows()-2,m.cols()-2))/pow(dy,2)
    +(m.block(1,2,m.rows()-2,m.cols()-2)+m.block(1,0,m.rows()-2,m.cols()-2))/pow(dx,2);
}


namespace {
  //returns a column of ones
  mesh2d::matrix onesCol(const int n) {return mesh2d::matrix::Ones(n,1);}

  //returns a row of ones
  mesh2d::matrix onesRow(const int n) {return mesh2d::matrix::Ones(1,n);}

  //copies first/last row/column from a to b.
  void copyBoundary(const mesh2d::matrix&a, mesh2d::matrix&b) {
    b.row(0)=a.row(0);
    b.col(0)=a.col(0);
    b.row(b.rows()-1)=a.row(a.rows()-1);
    b.col(b.cols()-1)=a.col(a.cols()-1);
    return;
  }

  //gives an RGB triple representing vl's position in the interval [min, max]
  //the return value will be black/white if vl is below/above the above interval.
  rgb_tuple rainbowScale(double _val, double min, double max) {
    const double val=5-(_val-min)/(max-min)*5;//scales from 0-6 if min<vl<max;
    constexpr float maxrgb=1;
    if (val<0) return rgb_tuple(0,0,0);//black if under range
    else if (val<1) return rgb_tuple(maxrgb,val*maxrgb,0);//red -> yellow
    else if (val<2) return rgb_tuple(maxrgb*(2-val),maxrgb,0);//yellow -> green
    else if (val<3) return rgb_tuple(0,maxrgb,maxrgb*(val-2));//green -> cyan
    else if (val<4) return rgb_tuple(0,maxrgb*(4-val),maxrgb);//cyan -> blue
    else if (val<=5) return rgb_tuple(maxrgb*(val-4),0,maxrgb);//blue -> magenta
    else return rgb_tuple(maxrgb,maxrgb,maxrgb);//white if val is outside usual range
  }
}


mesh2d::mesh2d(const int _nrow, const int _ncol, const double length, const double breadth):
  rho(1), mu(.01), CFL(.8), tolerance(.001),
  u_left(this->default_boundary), u_right(this->default_boundary), u_top(this->default_boundary), u_bottom(this->default_boundary),
  v_left(this->default_boundary), v_right(this->default_boundary), v_top(this->default_boundary), v_bottom(this->default_boundary),
  p_left(this->default_boundary), p_right(this->default_boundary), p_top(this->default_boundary), p_bottom(this->default_boundary)
{
  //rho, mu are density, viscosity resp.
  //Courant–Friedrichs–Lewy (CFL) criterion controls the time step;
  //the number CFL can be reduced if the solution diverges.
  double rho, mu, CFL;
  //the tolerance is used in a loop at solvePoisson
  double tolerance;
  nrow=abs(_nrow);ncol=abs(_ncol);
  //u, v, p are horizontal speed, vertical speed, pressure
  u=matrix(nrow+2,ncol+2);
  v=matrix(nrow+2,ncol+2);
  p=matrix(nrow+2,ncol+2);
  //u_star, v_star are solutions to the momentum equation ignoring pressure.
  u_star=matrix(nrow+2,ncol+2);
  v_star=matrix(nrow+2,ncol+2);
  //matrix derivatives of p, u_star, v_star
  p_xy=matrix(nrow,ncol);
  u_star_x=matrix(nrow,ncol);
  v_star_y=matrix(nrow,ncol);
  dx=length/(ncol-1);
  dy=breadth/(nrow-1);
  D=Differentiator(dx,dy);
  default_boundary=Boundary(Boundary::Dirichlet,0);
}

void mesh2d::setFluid(const double _rho, const double _mu, const double _CFL) {
  rho=_rho;
  mu=_mu;
  CFL=_CFL;
}

void mesh2d::setTolerance(const double _tolerance) {
  tolerance=_tolerance;
}

//finds an appropriate time step (dt) using the CFL criterion
void mesh2d::applyCFL () {
  const double denom=(u.array().maxCoeff()/dx+v.array().maxCoeff()/dy);
  dt=(denom==0)?CFL*(dx+dy):CFL/denom;
}

double mesh2d::getDt() {return dt;}

//Saves vertical velocity boundary conditions to the mesh2d structure
void mesh2d::setUBoundary(Boundary left, Boundary right, Boundary top, Boundary bottom) {
  u_left=left;u_right=right;u_top=top;u_bottom=bottom;
  resetUBoundary();
}

//Set boundary conditions for horizontal velocity
void mesh2d::resetUBoundary() {
  if (u_left.getType()==Boundary::Dirichlet) {
    u.col(0)=u_left.getValue()*onesCol(nrow+2);
  }
  else {
    u.col(0)=-u_left.getValue()*dx*onesCol(nrow+2)+u.col(1);
  }
  if (u_right.getType()==Boundary::Dirichlet) {
    u.col(ncol+1)=u_right.getValue()*onesCol(nrow+2);
  }
  else {
    u.col(ncol+1)=u_right.getValue()*dx*onesCol(nrow+2)+u.col(ncol);
  }
  if (u_top.getType()==Boundary::Dirichlet) {
    u.row(nrow+1)=2*u_top.getValue()*onesRow(ncol+2)-u.row(nrow);
  }
  else {
    u.row(nrow+1)=-u_top.getValue()*dy*onesRow(ncol+2)+u.row(nrow);
  }
  if (u_bottom.getType()==Boundary::Dirichlet) {
    u.row(0)=2*u_bottom.getValue()*onesRow(ncol+2)-u.row(1);
  }
  else {
    u.row(0)=u_bottom.getValue()*dy*onesRow(ncol+2)+u.row(1);
  }
}

//Saves vertical velocity boundary conditions to the mesh2d structure
void mesh2d::setVBoundary(Boundary left, Boundary right, Boundary top, Boundary bottom) {
  v_left=left;v_right=right;v_top=top;v_bottom=bottom;
  resetVBoundary();
}

//Set boundary conditions for vertical velocity
void mesh2d::resetVBoundary() {
  if (v_left.getType()==Boundary::Dirichlet) {
    v.col(0)=2*v_left.getValue()*onesCol(nrow+2)-v.col(1);
  }
  else {
    v.col(0)=-v_left.getValue()*dx*onesCol(nrow+2)+v.col(1);
  }
  if (v_right.getType()==Boundary::Dirichlet) {
    v.col(ncol+1)=2*v_right.getValue()*onesCol(nrow+2)-v.col(ncol);
  }
  else {
    v.col(ncol+1)=v_right.getValue()*dx*onesCol(nrow+2)+v.col(ncol);
  }
  if (v_top.getType()==Boundary::Dirichlet) {
    v.row(nrow+1)=v_top.getValue()*onesRow(ncol+2);
  }
  else {
    v.row(nrow+1)=-v_top.getValue()*dy*onesRow(ncol+2)+v.row(nrow);
  }
  if (v_bottom.getType()==Boundary::Dirichlet) {
    v.row(0)=v_bottom.getValue()*onesRow(ncol+2);
  }
  else {
    v.row(0)=v_bottom.getValue()*dy*onesRow(ncol+2)+v.row(1);
  }
}

//Saves pressure boundary conditions to the mesh2d structure for use in SolvePoisson(),
//then sets the boundary condition
void mesh2d::setPBoundary(Boundary left, Boundary right, Boundary top, Boundary bottom) {
  p_left=left;p_right=right;p_top=top;p_bottom=bottom;
  resetPBoundary();
}

//Set boundary conditions for pressure
void mesh2d::resetPBoundary() {
  if (p_left.getType()==Boundary::Dirichlet) {
    p.col(0)=p_left.getValue()*onesCol(nrow+2);
  }
  else {
    p.col(0)=-p_left.getValue()*dx*onesCol(nrow+2)+p.col(1);
  }
  if (p_right.getType()==Boundary::Dirichlet) {
    p.col(ncol+1)=p_right.getValue()*onesCol(nrow+2);
  }
  else {
    p.col(ncol+1)=p_right.getValue()*dx*onesCol(nrow+2)+p.col(ncol);
  }
  if (p_top.getType()==Boundary::Dirichlet) {
    p.row(nrow+1)=p_top.getValue()*onesRow(ncol+2);
  }
  else {
    p.row(nrow+1)=-p_top.getValue()*dy*onesRow(ncol+2)+p.row(nrow);
  }
  if (p_bottom.getType()==Boundary::Dirichlet) {
    p.row(0)=p_bottom.getValue()*onesRow(ncol+2);
  }
  else {
    p.row(0)=p_bottom.getValue()*dy*onesRow(ncol+2)+p.row(1);
  }
}

//performs an iteration of one time step
void mesh2d::doIteration() {
  applyCFL();
  //setting the boundary during each iteration is necessary if there are Neumann conditions.
  resetUBoundary();
  resetVBoundary();
  resetPBoundary();

  //updates u_star and v_star in the mesh2d structure
  //u_star, v_star are solutions to the momentum equation when pressure is ignored.
  auto u_face=(u.block(1,1,nrow,ncol)+u.block(1,2,nrow,ncol)
    +u.block(0,1,nrow,ncol)+u.block(0,2,nrow,ncol))/4;
  auto v_face=(v.block(1,1,nrow,ncol)+v.block(1,0,nrow,ncol)
    +v.block(2,1,nrow,ncol)+v.block(2,0,nrow,ncol))/4;
    copyBoundary(u,u_star);
    copyBoundary(v,v_star);
  u_star.block(1,1,nrow,ncol)=u.block(1,1,nrow,ncol)
    -dt*(u.block(1,1,nrow,ncol).cwiseProduct(D.x(u))+v_face.cwiseProduct(D.y(u)))
    +dt*(mu/rho)*(D.xx(u)+D.yy(u));
  v_star.block(1,1,nrow,ncol)=v.block(1,1,nrow,ncol)
    -dt*(u_face.cwiseProduct(D.x(v))+v.block(1,1,nrow,ncol).cwiseProduct(D.y(v)))
    +dt*(mu/rho)*(D.xx(v)+D.yy(v));

  //Solves the Poisson equation
  const double factor=1/(2/pow(dx,2)+2/pow(dy,2));
  u_star_x=D.x(u_star);
  v_star_y=D.y(v_star);
  double error=tolerance;
  for (int i=0;i<=500;++i) {
    if (tolerance>error) {break;}
    auto p_old=p;
    p_xy=D.xy(p);
    p.block(1,1,nrow,ncol)=p_xy*factor-
      //(rho*factor/dt)*(D_x(u_star)+D_y(v_star));
      (rho*factor/dt)*(u_star_x+v_star_y);
    error=(p-p_old).array().abs().maxCoeff();
    resetPBoundary();//necessary if there are Neumann conditions.
  }

  //Finally, solve the momentum equation
  u.block(1,1,nrow,ncol)=u_star.block(1,1,nrow,ncol)-(dt/rho)*D.x(p);
  v.block(1,1,nrow,ncol)=v_star.block(1,1,nrow,ncol)-(dt/rho)*D.y(p);
}

//writes to filename three tab-separated columns, giving the inner parts of the
//matrices p, u, v, flattened by concatenating successive rows.
void mesh2d::writeFile(const std::string& filename) {
  std::ofstream outfile;
  outfile.open (filename);
  for (int r=1;r<nrow+1;++r) {
    for (int c=1;c<ncol+1;++c) {
      outfile<<p(r,c)<<"\t"<<u(r,c)<<"\t"<<v(r,c)<<"\n";
    }
  }
  outfile.close();
  return;
}

/*********************************
This writes an image file to filename. Doubles min, max are used in calling
rainbow_scale above. The function
(1) creates an image with the same dimensions as the mesh2d object,
(2) colors each pixel by speed, based on rainbowScale
(3) places short lines (vectors) to represent the velocity direction (if the velocity is nonzero)
(4) places a white dot at the base of each vector
*********************************/
void mesh2d::writeImage(const std::string& filename, const double min, const double max) {
  //create a grid to place flow vectors
  std::vector<int> xpos, ypos;
  constexpr int width=7;//width/height of vectors
  constexpr int spacing=2;//spacing between vectors
  for (int i=width+spacing;i<nrow-spacing;i+=2*width+spacing) xpos.push_back(i);
  for (int i=width+spacing;i<ncol-spacing;i+=2*width+spacing) ypos.push_back(i);

  //find speed at eaceh point of the grid
  const matrix speed=(u.array().pow(2)+v.array().pow(2)).sqrt().matrix().block(1,1,nrow,ncol);

  //create image
  Magick::Image image(Magick::Geometry(nrow,ncol),Magick::Color(MaxRGB,MaxRGB,MaxRGB,0));

  //create rainbow-colored background to represent velocity at each grid point.
  rgb_tuple color;
  for (int i=0;i<nrow;++i) {
    for (int j=0;j<ncol;++j) {
      color=rainbowScale(speed.coeff(i,j),min,max);
      image.pixelColor(i,j,Magick::Color(color.r()*MaxRGB,color.g()*MaxRGB,color.b()*MaxRGB,MaxRGB));
    }
  }

  //add lines representing velocity vectors along a grid
  image.strokeColor("black");
  image.fillColor("black");
  image.strokeWidth(.5);
  std::list<Magick::Drawable> drawList;
  for (auto & x : xpos) {
    for (auto & y : ypos) {
      if (u.coeff(x+1,y+1)!=0 || v.coeff(x+1,y+1)!=0)
        drawList.push_back(Magick::DrawableLine(x,y,x+int(width*u.coeff(x+1,y+1)/speed.coeff(x,y)),
                                                y+int(width*v.coeff(x+1,y+1)/speed.coeff(x,y))));
    }
  }
  image.draw(drawList);

  //add a white dot at each point of the grid where we have placed a vector above
  for (auto & x : xpos) {
    for (auto & y : ypos) {
      image.pixelColor(x,y,Magick::Color(MaxRGB,MaxRGB,MaxRGB,MaxRGB));
    }
  }

  //write image to file
  try {
    image.write( filename );
  } catch( std::exception &error_ ) {
    std::cout << "Couldn't write to file "+filename << error_.what() << std::endl;
    return;
  }
}

double mesh2d::maxVel() {
  return (u.array().pow(2)+v.array().pow(2)).sqrt().block(1,1,nrow,ncol).maxCoeff();
}

double mesh2d::minVel() {
  return (u.array().pow(2)+v.array().pow(2)).sqrt().block(1,1,nrow,ncol).minCoeff();
}

mesh2d::matrix mesh2d::getU() {return u;}

mesh2d::matrix mesh2d::getV() {return v;}

mesh2d::matrix mesh2d::getP() {return p;}
