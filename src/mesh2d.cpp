#include "mesh2d.h"
#include <Magick++.h>
#include <iostream>
#include <fstream>

RGBTuple::RGBTuple() {}

RGBTuple::RGBTuple(const float _r, const float _g, const float _b) {
  red=_r;
  green=_g;
  blue=_b;
}

float RGBTuple::r() const {return red;}

float RGBTuple::g() const {return green;}

float RGBTuple::b() const {return blue;}


Boundary::Boundary(Condition _type, double _value) {
  value=_value;
  type=_type;
}

Boundary::Boundary(): type(Boundary::Condition::Dirichlet), value(0) {}

Boundary::Condition Boundary::getType() const {return type;}

double Boundary::getValue() const {return value;}


namespace {
  //returns a column of ones
  Mesh2D::matrix onesCol(const int n) {return Mesh2D::matrix::Ones(n,1);}

  //returns a row of ones
  Mesh2D::matrix onesRow(const int n) {return Mesh2D::matrix::Ones(1,n);}

  //copies first/last row/column from a to b.
  void copyBoundary(const Mesh2D::matrix&a, Mesh2D::matrix&b) {
    b.row(0)=a.row(0);
    b.col(0)=a.col(0);
    b.row(b.rows()-1)=a.row(a.rows()-1);
    b.col(b.cols()-1)=a.col(a.cols()-1);
    return;
  }

//Compute 1st derivative in x
  auto D_y(const Mesh2D::matrix& m, double dy) {
    return (m.block(2,1,m.rows()-2,m.cols()-2)-m.block(0,1,m.rows()-2,m.cols()-2))/(2*dy);
  }

  //Compute 1st derivative in x
  auto D_x(const Mesh2D::matrix& m, double dx) {
    return (m.block(1,2,m.rows()-2,m.cols()-2)-m.block(1,0,m.rows()-2,m.cols()-2))/(2*dx);
  }

  //Compute 2nd derivative in y
  auto D_yy(const Mesh2D::matrix& m, double dy) {
    return (m.block(2,1,m.rows()-2,m.cols()-2)-2*m.block(1,1,m.rows()-2,m.cols()-2)
      +m.block(0,1,m.rows()-2,m.cols()-2))/(pow(dy,2));
  }

  //Compute 2nd derivative in x
  auto D_xx(const Mesh2D::matrix& m, double dx) {
    return (m.block(1,2,m.rows()-2,m.cols()-2)-2*m.block(1,1,m.rows()-2,m.cols()-2)
      +m.block(1,0,m.rows()-2,m.cols()-2))/(pow(dx,2));
  }


  //Compute mixed partial derivative
  auto D_xy(const Mesh2D::matrix& m, double dx, double dy) {
      return (m.block(2,1,m.rows()-2,m.cols()-2)+m.block(0,1,m.rows()-2,m.cols()-2))/pow(dy,2)
        +(m.block(1,2,m.rows()-2,m.cols()-2)+m.block(1,0,m.rows()-2,m.cols()-2))/pow(dx,2);
  }

  //gives an RGB triple representing vl's position in the interval [min, max]
  //the return value will be black/white if vl is below/above the above interval.
  RGBTuple rainbowScale(double _val, double min, double max) {
    const double val=5-(_val-min)/(max-min)*5;//scales from 0-6 if min<vl<max;
    constexpr float maxrgb=1;
    if (val<0) return RGBTuple(0,0,0);//black if under range
    else if (val<1) return RGBTuple(maxrgb,val*maxrgb,0);//red -> yellow
    else if (val<2) return RGBTuple(maxrgb*(2-val),maxrgb,0);//yellow -> green
    else if (val<3) return RGBTuple(0,maxrgb,maxrgb*(val-2));//green -> cyan
    else if (val<4) return RGBTuple(0,maxrgb*(4-val),maxrgb);//cyan -> blue
    else if (val<=5) return RGBTuple(maxrgb*(val-4),0,maxrgb);//blue -> magenta
    else return RGBTuple(maxrgb,maxrgb,maxrgb);//white if val is outside usual range
  }
}


Mesh2D::Mesh2D(const int _nrow, const int _ncol, const double length, const double breadth):
  rho(1), mu(.01), CFL(.8), tolerance(.001),
  u_left(this->default_boundary), u_right(this->default_boundary),
    u_top(this->default_boundary), u_bottom(this->default_boundary),
  v_left(this->default_boundary), v_right(this->default_boundary),
    v_top(this->default_boundary), v_bottom(this->default_boundary),
  p_left(this->default_boundary), p_right(this->default_boundary),
    p_top(this->default_boundary), p_bottom(this->default_boundary)
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

  default_boundary=Boundary(Boundary::Condition::Dirichlet,0);
}

void Mesh2D::setFluid(const double _rho, const double _mu, const double _CFL) {
  rho=_rho;
  mu=_mu;
  CFL=_CFL;
}

void Mesh2D::setTolerance(const double _tolerance) {
  tolerance=_tolerance;
}

//finds an appropriate time step (dt) using the CFL criterion
void Mesh2D::applyCFL () {
  const double denom=(u.array().maxCoeff()/dx+v.array().maxCoeff()/dy);
  dt=(denom==0)?CFL*(dx+dy):CFL/denom;
}

double Mesh2D::getDt() const {return dt;}

//Saves vertical velocity boundary conditions to the Mesh2D structure
void Mesh2D::setUBoundary(const Boundary left, const Boundary right, const Boundary top, const Boundary bottom) {
  u_left=left;u_right=right;u_top=top;u_bottom=bottom;
  resetUBoundary();
}

//Set boundary conditions for horizontal velocity
void Mesh2D::resetUBoundary() {
  if (u_left.getType()==Boundary::Condition::Dirichlet) {
    u.col(0)=u_left.getValue()*onesCol(nrow+2);
  }
  else {
    u.col(0)=-u_left.getValue()*dx*onesCol(nrow+2)+u.col(1);
  }
  if (u_right.getType()==Boundary::Condition::Dirichlet) {
    u.col(ncol+1)=u_right.getValue()*onesCol(nrow+2);
  }
  else {
    u.col(ncol+1)=u_right.getValue()*dx*onesCol(nrow+2)+u.col(ncol);
  }
  if (u_top.getType()==Boundary::Condition::Dirichlet) {
    u.row(nrow+1)=2*u_top.getValue()*onesRow(ncol+2)-u.row(nrow);
  }
  else {
    u.row(nrow+1)=-u_top.getValue()*dy*onesRow(ncol+2)+u.row(nrow);
  }
  if (u_bottom.getType()==Boundary::Condition::Dirichlet) {
    u.row(0)=2*u_bottom.getValue()*onesRow(ncol+2)-u.row(1);
  }
  else {
    u.row(0)=u_bottom.getValue()*dy*onesRow(ncol+2)+u.row(1);
  }
}

//Saves vertical velocity boundary conditions to the Mesh2D structure
void Mesh2D::setVBoundary(const Boundary left, const Boundary right, const Boundary top, const Boundary bottom) {
  v_left=left;v_right=right;v_top=top;v_bottom=bottom;
  resetVBoundary();
}

//Set boundary conditions for vertical velocity
void Mesh2D::resetVBoundary() {
  if (v_left.getType()==Boundary::Condition::Dirichlet) {
    v.col(0)=2*v_left.getValue()*onesCol(nrow+2)-v.col(1);
  }
  else {
    v.col(0)=-v_left.getValue()*dx*onesCol(nrow+2)+v.col(1);
  }
  if (v_right.getType()==Boundary::Condition::Dirichlet) {
    v.col(ncol+1)=2*v_right.getValue()*onesCol(nrow+2)-v.col(ncol);
  }
  else {
    v.col(ncol+1)=v_right.getValue()*dx*onesCol(nrow+2)+v.col(ncol);
  }
  if (v_top.getType()==Boundary::Condition::Dirichlet) {
    v.row(nrow+1)=v_top.getValue()*onesRow(ncol+2);
  }
  else {
    v.row(nrow+1)=-v_top.getValue()*dy*onesRow(ncol+2)+v.row(nrow);
  }
  if (v_bottom.getType()==Boundary::Condition::Dirichlet) {
    v.row(0)=v_bottom.getValue()*onesRow(ncol+2);
  }
  else {
    v.row(0)=v_bottom.getValue()*dy*onesRow(ncol+2)+v.row(1);
  }
}

//Saves pressure boundary conditions to the Mesh2D structure for use in SolvePoisson(),
//then sets the boundary condition
void Mesh2D::setPBoundary(const Boundary left, const Boundary right, const Boundary top, const Boundary bottom) {
  p_left=left;p_right=right;p_top=top;p_bottom=bottom;
  resetPBoundary();
}

//Set boundary conditions for pressure
void Mesh2D::resetPBoundary() {
  if (p_left.getType()==Boundary::Condition::Dirichlet) {
    p.col(0)=p_left.getValue()*onesCol(nrow+2);
  }
  else {
    p.col(0)=-p_left.getValue()*dx*onesCol(nrow+2)+p.col(1);
  }
  if (p_right.getType()==Boundary::Condition::Dirichlet) {
    p.col(ncol+1)=p_right.getValue()*onesCol(nrow+2);
  }
  else {
    p.col(ncol+1)=p_right.getValue()*dx*onesCol(nrow+2)+p.col(ncol);
  }
  if (p_top.getType()==Boundary::Condition::Dirichlet) {
    p.row(nrow+1)=p_top.getValue()*onesRow(ncol+2);
  }
  else {
    p.row(nrow+1)=-p_top.getValue()*dy*onesRow(ncol+2)+p.row(nrow);
  }
  if (p_bottom.getType()==Boundary::Condition::Dirichlet) {
    p.row(0)=p_bottom.getValue()*onesRow(ncol+2);
  }
  else {
    p.row(0)=p_bottom.getValue()*dy*onesRow(ncol+2)+p.row(1);
  }
}

//performs an iteration of one time step
void Mesh2D::doIteration() {
  applyCFL();
  //setting the boundary during each iteration is necessary if there are Neumann conditions.
  resetUBoundary();
  resetVBoundary();
  resetPBoundary();

  //updates u_star and v_star in the Mesh2D structure
  //u_star, v_star are solutions to the momentum equation when pressure is ignored.
  auto u_face=(u.block(1,1,nrow,ncol)+u.block(1,2,nrow,ncol)
    +u.block(0,1,nrow,ncol)+u.block(0,2,nrow,ncol))/4;
  auto v_face=(v.block(1,1,nrow,ncol)+v.block(1,0,nrow,ncol)
    +v.block(2,1,nrow,ncol)+v.block(2,0,nrow,ncol))/4;
  copyBoundary(u,u_star);
  copyBoundary(v,v_star);
  u_star.block(1,1,nrow,ncol)=u.block(1,1,nrow,ncol)
    -dt*(u.block(1,1,nrow,ncol).cwiseProduct(D_x(u,dx))+v_face.cwiseProduct(D_y(u,dy)))
    +dt*(mu/rho)*(D_xx(u,dx)+D_yy(u,dy));
  v_star.block(1,1,nrow,ncol)=v.block(1,1,nrow,ncol)
    -dt*(u_face.cwiseProduct(D_x(v,dx))+v.block(1,1,nrow,ncol).cwiseProduct(D_y(v,dy)))
    +dt*(mu/rho)*(D_xx(v,dx)+D_yy(v,dy));

  //Solves the Poisson equation
  const double factor=1/(2/pow(dx,2)+2/pow(dy,2));
  u_star_x=D_x(u_star,dx);
  v_star_y=D_y(v_star,dy);
  double error=tolerance;
  for (int i=0;i<=500;++i) {
    if (tolerance>error) {break;}
    auto p_old=p;
    p_xy=D_xy(p,dx,dy);
    p.block(1,1,nrow,ncol)=p_xy*factor-
      (rho*factor/dt)*(u_star_x+v_star_y);
    error=(p-p_old).array().abs().maxCoeff();
    resetPBoundary();//necessary if there are Neumann conditions.
  }

  //Finally, solve the momentum equation
  u.block(1,1,nrow,ncol)=u_star.block(1,1,nrow,ncol)-(dt/rho)*D_x(p,dx);
  v.block(1,1,nrow,ncol)=v_star.block(1,1,nrow,ncol)-(dt/rho)*D_y(p,dy);
}

//writes to filename three tab-separated columns, giving the inner parts of the
//matrices p, u, v, flattened by concatenating successive rows.
void Mesh2D::writeFile(const std::string& filename) const {
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
rainbowScale above. The function
(1) creates an image with the same dimensions as the Mesh2D object,
(2) colors each pixel by speed, based on rainbowScale
(3) places short lines (vectors) to represent the velocity direction (if the velocity is nonzero)
(4) places a white dot at the base of each vector
*********************************/
void Mesh2D::writeImage(const std::string& filename, const double min, const double max) const {
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
  RGBTuple color;
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

const Mesh2D::matrix& Mesh2D::getU() const {return u;}

const Mesh2D::matrix& Mesh2D::getV() const {return v;}

const Mesh2D::matrix& Mesh2D::getP() const {return p;}
