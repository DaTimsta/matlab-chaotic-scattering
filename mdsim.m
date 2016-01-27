% molecular dynamics
% mdsim(lattice,n,dt,xp,yp,theta,v0)
% create lattice with create_lattice.m
% n is number of steps
% dt is time step
% xp,yp is starting position for particle
% v0,theta is starting velocity for particle
% outputs the phase space r,v
% 
% symbolic version
% takes the symbolic potential U depending on symbolic variables
% x and y

function [r,v] = mdsim(n,dt,xp,yp,theta,v0,U,x,y)

% create initial particle
r = zeros(n+1,2);
r(1,1) = xp;
r(1,2) = yp;

%with initial velocity
v = zeros(n+1,2);
v(1,1) = v0*cos(theta);
v(1,2) = v0*sin(theta);

% initialize r_prev
r_prev = r(1,:) - v(1,:).*dt;

forceX_mlab = matlabFunction(-diff(U,x));
forceY_mlab = matlabFunction(-diff(U,y));

for i=1:n %time
    % force is the minus the gradient of the potential
    % so we integrate directly
    force = [forceX_mlab(r(i,1),r(i,2)) forceY_mlab(r(i,1),r(i,2))];
    r_new = 2.*r(i,:) - r_prev + force.*(dt.^2);
    v(i+1,:) = (r_new - r_prev)/(2*dt);
    r_prev = r(i,:);
    r(i+1,:) = r_new;
end