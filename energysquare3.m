function energy = energysquare3(Sx,Sy,Sz,Jx,Jy,Dx,Dy,K)
%energy of a rectangular with periodic boundary conditions
%Exchange interactions
%Jx,Jy > 0 are nn exchange constants in the x,y directions
%The negative sign for AFM is taken into account later

%DM interactions
%x-direction
% Dx*[nz d/dx ny - ny d/dx nz]
% 0.5*Dx*[Sz(n)*(Sy(n+x) - Sy(n-x)) - Sy(n)*(Sz(n+x) - Sz(n-x))]
%pbc:  Dx*[Sz(n)*Sy(n+x) - Sz(n+x)*Sy(n)]

%y-direction
% Dy*[nx d/dy nz - nz d/dy nx]
% 0.5*Dy*[Sx(n)*(Sz(n+y) - Sz(n-y)) - Sz(n)*(Sx(n+y) - Sx(n-y))]
%pbc:  Dy*[Sx(n) Sz(n+y) - Sx(n+y)*Sz(n)]
%

%K is anisotropy (E = K/2Sz^2)
[Ny,Nx] = size(Sx);

%Useful ranges

%Ix = 1:Nx; 
%Iy = 1:Ny; 
IPx = [2:Nx 1]; %shifted index for nearest neighbours
%IPPx = [3:Nx 1 2];
%IMx = [Nx 1:Nx-1]; 
%IMMx = [Nx-1 Nx 1:Nx-2];

IPy = [2:Ny 1];
%IPPy = [3:Ny 1 2];
%IMy = [Ny 1:Ny-1]; 
%IMMy = [Ny-1 Ny 1:Ny-2];

%nn Heisenberg exchange
E1 =   -Sx.*(Jx*Sx(:,IPx)+Jy*Sx(IPy,:))...
       -Sy.*(Jx*Sy(:,IPx)+Jy*Sy(IPy,:))...
       -Sz.*(Jx*Sz(:,IPx)+Jy*Sz(IPy,:));
e1 = sum(sum(E1));

%DM interaction
E2 = Dx*(Sz.*Sy(:,IPx) - Sz(:,IPx).*Sy)...
    +Dy*(Sx.*Sz(IPy,:) - Sx(IPy,:).*Sz);
e2 = sum(sum(E2));

%anisotropy
e3 = -0.5*K*sum(sum(Sz.^2-1.0));

%E at inf
EFM = -(Jx+Jy);

%Final energy
energy = e1 + e2 + e3 - Nx*Ny*EFM;

end
