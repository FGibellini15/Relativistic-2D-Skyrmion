%Function to evaluate value of Neel vector at given point,
%according to the analytical profile discussed in the paper

function [nx,ny,nz] = RelativisticProfile_tweaked(x,y,X,Y,beta,lx,ly,R,w,A)
%The function needs as input a lattice (x,y), the coordinates of the lattice
%center, the speed of the skyrmion, the additional contraction and the
%parameters defining the size of the skyrmion

%Relativistic coefficient
gamma = 1/sqrt(1-beta^2);

%Coordinates, with proper contraction
deltax = gamma*lx*(x-X); %Distance from center
deltay = ly*(y-Y);
r = sqrt(deltax.^2 + deltay.^2);

%Defining angles
theta = 2*atan(sinh(R/w)./sinh(r/w));

phi = atan2(deltay,deltax);
phi = phi + pi/2 + A*sin(2*phi); %slight tweaking for better fitting

%Defining vector components
nx = cos(phi) .* sin(theta);
ny = sin(phi) .* sin(theta);
nz = cos(theta);

end

