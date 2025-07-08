function [HeffX,HeffY,HeffZ]=EffFieldsquare3(Sx,Sy,Sz,Jx,Jy,Dx,Dy,K)
%effective field for a square lattice (see energysquare3.m)
%periodic boundary conditions
%Jx,Jy > 0 (FM) are nn exchange constants in the x,y directions

%DM interactions

%x-direction
% Dx*[nz d/dx ny - ny d/dx nz]
% 0.5*Dx*[Sz(n)*(Sy(n+x) - Sy(n-x)) - Sy(n)*(Sz(n+x) - Sz(n-x))]
%pbc:  Dx*[Sz(n)*Sy(n+x) - Sz(n+x)*Sy(n)]
%Hy = -Dx*(Sz(n-x) - Sz(n+x))
%Hz = -Dx*(Sy(n+x) - Sy(n-x))

%y-direction
% Dy*[nx d/dy nz - nz d/dy nx]
% 0.5*Dy*[Sx(n)*(Sz(n+y) - Sz(n-y)) - Sz(n)*(Sx(n+y) - Sx(n-y))]
%pbc:  Dy*[Sx(n) Sz(n+y) - Sx(n+y)*Sz(n)]
%Hz = -Dy*(Sx(n-y) - Sx(n+y))
%Hx = -Dy*(Sz(n+y) - Sz(n-y))

%K is anisotropy (E = K/2Sz^2)
[Ny,Nx] = size(Sx);
%Ix = 1:Nx; 
%Iy = 1:Ny; 
IPx = [2:Nx 1];
%IPPx = [3:Nx 1 2];
IMx = [Nx 1:Nx-1]; 
%IMMx = [Nx-1 Nx 1:Nx-2];
%
IPy = [2:Ny 1];
%IPPy = [3:Ny 1 2];
IMy = [Ny 1:Ny-1]; 
%IMMy = [Ny-1 Ny 1:Ny-2];
%effective field
HeffX = Jx*(Sx(:,IPx)+Sx(:,IMx))+Jy*(Sx(IPy,:)+Sx(IMy,:))...
       -Dy*(Sz(IPy,:) - Sz(IMy,:)); 
%
HeffY = Jx*(Sy(:,IPx)+Sy(:,IMx))+Jy*(Sy(IPy,:)+Sy(IMy,:))...
       -Dx*(Sz(:,IMx) - Sz(:,IPx));  
%Hy = -Dx*(Sz(n-x) - Sz(n+x))

HeffZ = Jx*(Sz(:,IPx)+Sz(:,IMx))+Jy*(Sz(IPy,:)+Sz(IMy,:))...
       -Dx*(Sy(:,IPx) - Sy(:,IMx))-Dy*(Sx(IMy,:)-Sx(IPy,:))...
       +K*Sz;
%Hz = -Dx*(Sy(n+x) - Sy(n-x))

end