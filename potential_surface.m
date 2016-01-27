% plot the potential surface given a lattice and a symbolic potential
clear all
close all

% lattice parameters
N = [2 2];
d = 3;
% center the lattice
x0 = -d*(N(1)-1)*0.5;
y0 = -d*(N(2)-1)*0.5;

% create r^6 wall?
wall = 0;

% define our potential symbolically
syms x y lx ly
P = 4*(1/((x-lx)^2 + (y-ly)^2)^3)*((1/((x-lx)^2 + (y-ly)^2)^3)-1);

% create lattice
lattice = create_lattice(N(1),N(2),d,x0,y0);

% get the specific potential U
[U,x,y] = sym_potential(lattice,P,x,y,lx,ly);

% % get the jacobian
% gradU = [diff(U,x);diff(U,y)];
% JU = [diff(gradU(1),x) diff(gradU(1),y); diff(gradU(2),x) diff(gradU(2),y)];
% 
% % find jacobian at 0,0
% JU_fp = subs(subs(JU,x,-1.188123431174824),y,0);
% 
% disp(eig(JU_fp));

spacing = 1e-1; %precision of plot

U_mlab = matlabFunction(U);

X = (x0-d):spacing:(lattice(end,1)+d);
Y = (y0-d):spacing:(lattice(end,2)+d);
[Xg,Yg] = meshgrid(X,Y);

Z = U_mlab(Xg,Yg);

% constraint on the  graph
Z(Z>1) = 1;

figure
surf(Xg,Yg,Z);

figure
contour(Xg,Yg,Z,30);