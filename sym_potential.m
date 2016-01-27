% [U,x,y] = sym_potential(lattice,P)
% given a lattice using create_lattice.m
% and a potential P (depending on x,y,lx,ly)
% where lx and ly is distance to the lattice points
% create a symbolic lattice U depending on x and y 
% in cartesian coordinates
function [U,x,y] = sym_potential(lattice,P,x,y,lx,ly)
U = 0; %init

% create the total potential as a sum of its parts
for i=1:length(lattice)
    dlx = lattice(i,1);
    dly = lattice(i,2);
    U = U + subs(subs(P,ly,dly),lx,dlx);
end