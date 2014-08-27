function [dudx,dudy,dvdx,dvdy,vor,div,vbyf] = strainrateest(u,v,x,y,Rmax,Rmin,Nmax);
#function [dudx,dudy,dvdx,dvdy,vor,div,vbyf] = strainrateest(u,v,x,y,Rmax,Rmin,Nmax);
# given a set of velocities and locations calculate the derivatives
# using a least squares aproximation
# u = zonal velocity (m/s)
# v = meridional velocity (m/s)
# x = decimal longitude in degrees
# y = decimal latitude in degrees
# Rmax = point used in the derivative must be closer that this (default 12 000 m)
# Rmin = point used in the derivative must be further that this (default 1 000 m)
# N = maximum number of points to consider (default is 12)
# uses the Gibbs sea water package
 addpath /Users/mhoecker/Documents/TEOS-10
 addpath /Users/mhoecker/Documents/TEOS-10/library
# Set default value for Rmax
 if(nargin<5)
  Rmax = 12000;
 end%if
# Set default value for Rmin
 if(nargin<6)
  Rmin = 1000;
 end%if
# Set default value for Nmax
 if(nargin<7)
  Nmax = 12;
 end%if
# Initialize derivitive lists
 dudx = zeros(size(u));
 dudy = zeros(size(u));
 dvdx = zeros(size(v));
 dvdy = zeros(size(v));
#cycle over points
 for j=1:length(u)
# Define the reference point
  x0 = x(j);
  y0 = y(j);
  u0 = u(j);
  v0 = v(j);
  dist = zeros(size(x));
  for i=1:length(x)
   dist(i) = gsw_distance([x(i),x0],[y(i),y0]);
  end%for
  idxnear = find((dist<Rmax).*(dist>Rmin).*(isnan(dist+u+v+x+y)==0));
  N = length(idxnear);
# Ensure there are enough points and the reference is not NaNs
  if((N>3).*(isnan(x0+y0+u0+v0)==0))
   dx = zeros(2*N,4);
   du = zeros(2*N,1);
# These matricies form the equation
#
# dx*der = du
#
#
# where
#
# dx =
# dx1 dy1 000 000
# 000 000 dx1 dy1
# dx2 dy2 000 000
# 000 000 dx2 dy2
# .. .. .. ..
#
# du matrix
# du1
# dv1
# du2
# dv2
# .
# .
#
# der=
# du/dx
# du/dy
# dv/dx
# dv/dy
#
   for i=1:N
# x-x0
    dx(2*i-1,1) = gsw_distance([x(idxnear(i)),x0],[y0,y0])*sign(x(idxnear(i))-x0);
# y-y0
    dx(2*i-1,2) = gsw_distance([x0,x0],[y(idxnear(i)),y0])*sign(y(idxnear(i))-y0);
# u-u0
    du(2*i-1)   = u(idxnear(i))-u0;
# x-x0
    dx(2*i,3)   = gsw_distance([x(idxnear(i)),x0],[y0,y0])*sign(x(idxnear(i))-x0);
# y-y0
    dx(2*i,4)   = gsw_distance([x0,x0],[y(idxnear(i)),y0])*sign(y(idxnear(i))-y0);
# v-v0
    du(2*i)     = v(idxnear(i))-v0;
   end%for
   der =	ols(du,dx);
   dudx(j) = der(1);
   dudy(j) = der(2);
   dvdx(j) = der(3);
   dvdy(j) = der(4);
  else
# set derivatives to NaN if the point is a NaN or there are not enought nearby points
   dudx(j) = NaN;
   dudy(j) = NaN;
   dvdx(j) = NaN;
   dvdy(j) = NaN;
  end%if
end%for
# divergence = du/dx+dv/dy
div = dudx+dvdy;
# vorticity = dv/dx-du/dy
vor = dvdx-dudy;
# Scaled vorticity
vbyf = vor./gsw_f(y);
# THE END :)
end%function