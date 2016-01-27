% simulations using symbolic potential
clear all
close all

% steps (total time n*dt)
n = 3e5;
dt = 0.001;

% lattice parameters
N = [2 2];
d = 1.3;
x0 = 0;
y0 = 0;

% particle parameters
xp = d-0.01*d;
yp = d+0.02*d;
theta = pi/2-0.15;
v0 = 0;

% define our potential symbolically
syms x y lx ly
P = 4*(1/((x-lx)^2 + (y-ly)^2)^3)*((1/((x-lx)^2 + (y-ly)^2)^3)-1);

% lets vary the initial conditions
figure
hold on

ic_num = 1; % number of initial conditions
colors = colormap(hsv(ic_num));
lattice = create_lattice(N(1),N(2),2*d,x0,y0);

% get the specific potential U
[U,x,y] = sym_potential(lattice,P,x,y,lx,ly);

for i=1:ic_num
    [r,v] = mdsim(n,dt,xp,yp+0.001*i*d,theta,v0,U,x,y);
    plot(r(:,1),r(:,2),'color',colors(i,:),'LineWidth',2.5); %trajectory 
end

plot(lattice(:,1),lattice(:,2),'k*','MarkerSize',14,'LineWidth',2);
plot(xp,yp,'ko','MarkerSize',14,'LineWidth',2);
xlabel('x')
ylabel('y')
xlim([-d 3*d]);
ylim([-d 3*d]);
title('d=1.3');

U_mlab = matlabFunction(U);
% energy
disp(1/2.*(v(1,1).^2 + v(1,2).^2) + U_mlab(r(1,1),r(1,2)));

% power spectrum
N = length(r)-1; %% number of points
T = dt*n; %% define time of interval, 3.4 seconds
t = (0:N-1)/N; %% define time
t = t*T; %% define time in seconds
f = r(1:end-1); %%define function, 10 Hz sine wave
p = abs(fft(f))/(N/2); %% absolute value of the fft
p = p(1:N/2).^2; %% take the power of positve freq. half
freq = (0:N/2-1)/T; %% find the corresponding frequency in Hz
figure
semilogx(freq,p,'.-'); %% plot on semilog scale