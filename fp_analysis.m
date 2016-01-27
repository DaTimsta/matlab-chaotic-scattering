% fixed point analysis along the symmetric axes.

% plot the potential surface given a lattice and a symbolic potential
clear all
%close all

% bifurcations at d=1.05,1.2,1.3
d = 1.05;
title_string = 'Fixed points of model with d=1.3';
r_values = d/2:1e-2:3*d/2;
r_vector = zeros(1,length(r_values));

syms y
U = 8*(1/(d^2+(y-d)^2)^3 - 1/(d^2+(y-d)^2)^6) + 8*(1/(d^2+(y+d)^2)^3 - 1/(d^2+(y+d)^2)^6);
dU = diff(U,y);
dU_mlab = matlabFunction(dU);
d2U = diff(dU,y);
d2U_mlab = matlabFunction(d2U);

% newton-rhapson init, with guesses
ep = 1e-10;
n_max = 100;

for j=1:length(r_values);

    r = r_values(j);

    for i=1:n_max
        r_new = r - dU_mlab(r)./d2U_mlab(r);
        error = abs(dU_mlab(r_new));
        if error < ep
            %disp(i)
            %disp(error)
            %disp(r_new);
            %disp(d);
            r_vector(j) = r_new;
            break
        end

        r = r_new;
    end
end

% find the unique values to the precision defined by newtons method
% by using round() and unique(). Quite ugly code, but it seems to work.
fixed_points = ep*unique(round(r_vector(r_vector < 2*d & r_vector > 0).*(1/ep)));
disp(fixed_points);

% figuring out the nature of the fixed points

% having found the symmetric fixed points, 
% we will now draw them on a contour plot
% this code is a bit messy because of overlap from code used in other
% places

% lattice parameters
N = [2 2];
d = 2*d;
% center the lattice
x0 = -d*(N(1)-1)*0.5;
y0 = -d*(N(2)-1)*0.5;

% define our potential symbolically
syms x y lx ly
P = 4*(1/((x-lx)^2 + (y-ly)^2)^3)*((1/((x-lx)^2 + (y-ly)^2)^3)-1);

% create lattice
lattice = create_lattice(N(1),N(2),d,x0,y0);

% get the specific potential V
[V,x,y] = sym_potential(lattice,P,x,y,lx,ly);
V_mlab = matlabFunction(V);
L11 = diff(diff(V,x),x); L11_mlab = matlabFunction(L11);
L12 = diff(diff(V,x),y); L12_mlab = matlabFunction(L12);
L21 = diff(diff(V,y),x); L21_mlab = matlabFunction(L21);
L22 = diff(diff(V,y),y); L22_mlab = matlabFunction(L22);

spacing = 1e-2; %precision of plot

X = (x0-d):spacing:(lattice(end,1)+d);
Y = (y0-d):spacing:(lattice(end,2)+d);
[Xg,Yg] = meshgrid(X,Y);

Z = V_mlab(Xg,Yg);

% constraint on the  graph
Z(Z>0.5) = 0.5;

figure
hold on
contour(Xg,Yg,Z,15);
xlabel('x')
ylabel('y')
%title(title_string);

% plot the fixed points
for k=1:length(fixed_points)
    L = [L11_mlab(0,fixed_points(k)) L12_mlab(0,fixed_points(k));
        L21_mlab(0,fixed_points(k)) L22_mlab(0,fixed_points(k))];
    L_ev = eig(L); % Eigenvalues of L
    disp(L_ev);
    if sign(L_ev(1)) == sign(L_ev(2))        
        % center point, plot with a o
        plot(0,fixed_points(k),'ko','MarkerSize',10,'LineWidth',2.5);
        plot(0,-fixed_points(k),'ko','MarkerSize',10,'LineWidth',2.5);
        plot(fixed_points(k),0,'ko','MarkerSize',10,'LineWidth',2.5);
        plot(-fixed_points(k),0,'ko','MarkerSize',10,'LineWidth',2.5);
    else
        % saddle point, plot with an x
        plot(0,fixed_points(k),'kx','MarkerSize',12,'LineWidth',2.5);
        plot(0,-fixed_points(k),'kx','MarkerSize',12,'LineWidth',2.5);
        plot(fixed_points(k),0,'kx','MarkerSize',12,'LineWidth',2.5);
        plot(-fixed_points(k),0,'kx','MarkerSize',12,'LineWidth',2.5);
    end
end;