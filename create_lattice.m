% lattice = create_lattice(N1,N2,d,x0,y0)
% creates lattice with (N1)x(N2) points, distance d between them
% lower-right lattice point will start at x0,y0
function lattice = create_lattice(N1,N2,d,x0,y0)

lattice_size = N1*N2;
lattice = zeros(lattice_size,2);

for i=1:N1
    for j=1:N2
        lattice_num = (i-1)*N2 + j;
        lattice(lattice_num,1) = x0 + (i-1)*d;
        lattice(lattice_num,2) = y0 + (j-1)*d;
    end
end
