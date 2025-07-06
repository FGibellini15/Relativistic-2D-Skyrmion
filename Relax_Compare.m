close all

%Loading texture
load skyrmion5_61.mat

%Obtaining useful variables
[Sx,Sy,Sz] = normalizeSpins(Sx,Sy,Sz);
[Nx,Ny] = size(Sx);
[x,y] = meshgrid(0:1:Nx-1,0:1:Ny-1);
X = floor((Nx-1)/2);
Y = floor((Ny-1)/2);

%Defining desired parameters
J = 1;
K = 0.05;
beta = 0.58;
rho = 0.98;

%Basic and Relativistic coefficients
gamma = 1./sqrt(1-beta.^2);
D = rho*2*sqrt(K*J)/pi;
[lx,ly] = ContractionCompute_raw(J,K,rho,beta);
%lx = 1;
%ly = 1;

%Interaction Coefficients along directions
Jx = J./(gamma).*lx./ly;         %interaction along x
Jy = J./(gamma).*ly./lx;        %interaction along y
Dx = D./ly; %0.17;   %DM strength along x
Dy = D./gamma.*(1./lx);        %DM strength along y
K =  K./(gamma).*1./(lx.*ly);      %K > 0 is an easy-axis anisotropy

%Relaxing the numerical skyrmion and computing its energy
[Sx,Sy,Sz] = Relax(Sx,Sy,Sz,Jx,Jy,Dx,Dy,K,8000,2000,0.1);
energysquare3(Sx,Sy,Sz,Jx,Jy,Dx,Dy,K);

%Fitting the analytical profile
[R,w,A,fmin] = MinimizeDeltaE(Sx,Sy,Sz,J,K,beta,rho)
[nx,ny,nz] = RelativisticProfile_tweaked(x,y,X,Y,beta,lx,ly,R,w,A);

%Producing quiver plots of the skyrmions
Texture_plot(Nx,Ny,Sx,Sy,Sz);
title(['\beta = ' num2str(beta) ' \rho = ' num2str(rho)], FontSize=15 )
%title('Numeric')

%Texture_plot(Nx,Ny,nx,ny,nz);
%title('Fitted')