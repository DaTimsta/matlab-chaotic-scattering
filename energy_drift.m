% energy drift test
clear all
close all

% steps (total time n*dt)
n = 2e4;
dt = 0.001;

% lattice parameters
N = [3 3];
d = 1.2;
x0 = 2.5;
y0 = 0;

% particle parameters
xp = 0;
yp = d+0.14021;
theta = 0;
E = 2; % starting kinetic energy
v0 = sqrt(2*E);

% define our potential symbolically
syms x y lx ly
P = 4*(1/((x-lx)^2 + (y-ly)^2)^3)*((1/((x-lx)^2 + (y-ly)^2)^3)-1);
lattice = create_lattice(N(1),N(2),2*d,x0,y0);
[U,x,y] = sym_potential(lattice,P,x,y,lx,ly);

% simulate
[r,v] = mdsim(n,dt,xp,yp,theta,v0,U,x,y);

% plot trajectory and lattice
figure
hold on
plot(r(:,1),r(:,2),'LineWidth',3);
plot(lattice(:,1),lattice(:,2),'k*','MarkerSize',14,'LineWidth',2);
xlabel('x');
ylabel('y');
title('Trajectory associated with energy drift');


% energy considerations
U_mlab = matlabFunction(U);
etot = 1/2.*(v(:,1).^2 + v(:,2).^2) + U_mlab(r(:,1),r(:,2));

figure
plot(0:dt:n*dt,etot,'LineWidth',2.2);
ylabel('Total energy')
xlabel('Time')
title('Energy drift')
disp(std(etot));
