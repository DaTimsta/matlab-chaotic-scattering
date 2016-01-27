%scattering function analysis.
clear all
close all

% steps (total time n*dt)
n = 5e4;
dt = 0.001;

% lattice parameters
N = [4 4];
d = 1.3;
x0 = 2.5;
y0 = 0;

% particle parameters
xp = 0;
yp = 0;
theta = 0;
v0 = 2;

b_num = 10000; %number of times we change the impact parameter
precision =6*d/b_num; %the resulting precision
b_values = linspace(0,6*d,b_num);
phi = zeros(1,length(b_num));

% define our potential symbolically
syms x y lx ly
P = 4*(1/((x-lx)^2 + (y-ly)^2)^3)*((1/((x-lx)^2 + (y-ly)^2)^3)-1);

lattice = create_lattice(N(1),N(2),2*d,x0,y0);

% get the specific potential U
[U,x,y] = sym_potential(lattice,P,x,y,lx,ly);
r_end = zeros(b_num,2);
%track progress, simulations take a while
waitbar_handle = waitbar(0,'Initializing waitbar...');

for i=1:b_num
    tic;
    [r,v] = mdsim(n,dt,xp,precision*i,theta,v0,U,x,y);
    r_end(i,:) = [r(end,1);r(end,2)];
    phi(i) = sign(v(end,2))*acos(v(end,1)/sqrt(v(end,1)^2 + v(end,2)^2));
    waitbar(i/b_num,waitbar_handle,sprintf('%d out of %d done...',i,b_num));
    toc;
end

close(waitbar_handle);

%% Plots. Requires the following data: phi, lattice, r_end, b_values
close all
 
% plot endpoints and lattice as a check
figure
hold on
plot(lattice(:,1),lattice(:,2),'k*','MarkerSize',10);
plot(r_end(:,1),r_end(:,2),'bo');

% plot the scattering function
figure
plot(b_values,phi,'.');
title('Scattering function with d=1.3');
xlabel('b');
ylabel('\phi');
axis tight

%scattering cross section
resolution = [0.018;0.009];

for j=1:length(resolution)
    dcs_size = floor(2*pi/resolution(j));
    dcs = zeros(0,dcs_size);
    for k=1:dcs_size
        dcs(k) = length(phi(phi > resolution(j)*k - pi & phi < resolution(j)*(k+1) - pi));
    end
    
    figure
    plot(linspace(-pi,pi,dcs_size),dcs,'-','LineWidth',2)
    xlabel('\phi');
    ylabel('diff. cross section');
    axis tight
end

%%