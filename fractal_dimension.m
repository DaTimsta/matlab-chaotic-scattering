% fractal_dimension
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
b_domain = [3.5;4.3]; % look at the fractal part of the graph

% particle parameters
xp = 0;
yp = 0;
theta = 0;
v0 = 2;

% define our potential symbolically
syms x y lx ly
P = 4*(1/((x-lx)^2 + (y-ly)^2)^3)*((1/((x-lx)^2 + (y-ly)^2)^3)-1);
lattice = create_lattice(N(1),N(2),2*d,x0,y0);

% get the specific potential U
[U,x,y] = sym_potential(lattice,P,x,y,lx,ly);

epsilon = [5e-8;1e-7;5e-7;1e-6;5e-6;1e-5;5e-5;1e-4;5e-4;1e-3;5e-3];
b_num = 500; %number of times we change the impact parameter per epsilon value
f_epsilon = zeros(1,length(epsilon));

%track progress, simulations take a while
%waitbar_handle = waitbar(0,'Initializing waitbar...');

for i=1:length(epsilon)
    count = 0;
    for j=1:b_num
        tic;
        b1 = rand()*diff(b_domain) + b_domain(1);
        if randi(2) == 1            
            b2 = b1 + epsilon(i);
        else
            b2 = b1 - epsilon(i);
        end
        [r,v] = mdsim(n,dt,xp,b1,theta,v0,U,x,y);
        s1 = sign(v(end,2));
        [r,v] = mdsim(n,dt,xp,b2,theta,v0,U,x,y);
        s2 = sign(v(end,2));
        
        if s1 ~= s2
            count = count + 1;
        end        
        %waitbar(j/b_num,waitbar_handle,sprintf('%d out of %d done...',i,length(epsilon)));
        toc;
    end
    f_epsilon(i) = count/b_num;
end

%%

close all
%find fractal dimension
fractal_fit = polyfit(log(epsilon),log(f_epsilon'),1);
alpha = fractal_fit(1);
beta = fractal_fit(2);

figure
loglog(epsilon,f_epsilon,'rd','MarkerSize',10,'LineWidth',2);
hold on
loglog(epsilon,exp(beta + log(epsilon)*alpha),'b-','LineWidth',2.5);
%plot(log(epsilon),log(f_epsilon),'rd','MarkerSize',10,'LineWidth',2);
%plot(log(epsilon),beta + log(epsilon)*alpha,'b-','LineWidth',3);
xlabel('\epsilon')
ylabel('f(\epsilon)');

disp(alpha)