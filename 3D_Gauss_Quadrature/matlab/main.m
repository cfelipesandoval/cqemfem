clear; clc

% 1 tetra
f_kernel_tetra = @(x,y,z) x+y+z;
nodes_tetra1 = rand(4,3);
quadTetra(f_kernel_tetra, nodes_tetra1)

% 2 tetra
%k = 1;
%g = @(xm,ym,zm,xn,yn,zn) (exp(1i*k*sqrt((xm-xn).^2+(ym-yn).^2+(zm-zn).^2)) /4/pi./ sqrt((xm-xn).^2+(ym-yn).^2+(zm-zn).^2) );
%nodes_tetra1 = rand(4,3);
%nodes_tetra2 = rand(4,3) + 10;
%quadTetraTetra(g, nodes_tetra1, nodes_tetra2)