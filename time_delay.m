% calculation of mean exit time
clear all
close all

% steps (total time n*dt)
n = 5e4;
dt = 0.001;

% lattice parameters (remember to center)
N = [4 4];
d = 1.3;
x0 = -3*d;
y0 = -3*d;
lattice = create_lattice(N(1),N(2),2*d,x0,y0);

% define our potential symbolically
syms x y lx ly
P = 4*(1/((x-lx)^2 + (y-ly)^2)^3)*((1/((x-lx)^2 + (y-ly)^2)^3)-1);
[U,x,y] = sym_potential(lattice,P,x,y,lx,ly);

% chose how we sprinkle the initial conditions about the origin.
v0 = 2; % held constant
num_ic = 40;
yp = 0;
xp = linspace(-0.5,-0.9,num_ic);
theta = linspace(1.4,1.8,num_ic);

% size of the region R (radial)
R_size = 2.5 + d/sqrt(2);
% initialize vector that contains the times particle leaves
left_R = zeros(num_ic,num_ic);

%track progress, simulations take a while
waitbar_handle = waitbar(0,'Initializing waitbar...');

%figure
%hold on

for i=1:num_ic
    for j=1:num_ic
        tic;
        [r,v] = mdsim(n,dt,xp(i),yp,theta(j),v0,U,x,y);
        %plot(r(:,1),r(:,2))
        % figure out when we leave the scattering region
        for k=1:length(r)
            if sqrt(r(k,1)^2 + r(k,2)^2) > R_size
                left_R(i,j) = k;
                break;
            end
        end
        
        waitbar(j/num_ic,waitbar_handle,sprintf('%d out of %d done...',i,num_ic));
        toc;
    end
end

close(waitbar_handle);

left_R = left_R ( : );
N_T = zeros(1,n+1);
N_T(1) = num_ic^2;

for l=1:n
    N_T(l+1) = N_T(l) - length(left_R(left_R == l));
end
