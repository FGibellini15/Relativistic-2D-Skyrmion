function deltaE = energy_diff(Sx,Sy,Sz,J,K,beta,rho,lx,ly,R,w,A)
    %This function computes the absolute value of the energy difference 
    %between the given texture and the analytical profile evaluated at
    %given speed and rho

    %Useful quantities of the lattice
    [Nx,Ny] = size(Sx);
    [x,y] = meshgrid(0:1:Nx-1,0:1:Ny-1);
    X = floor((Nx-1)/2); %Coordinates of the center
    Y = floor((Ny-1)/2);

    %Necessary coefficients
    gamma = 1./sqrt(1-beta.^2); %relativistic
    D = rho*2*sqrt(K*J)/pi;     %DM interaction strength

    %Interaction Coefficients
    Jx = J./(gamma).*lx./ly;         %interaction along x
    Jy = J./(gamma).*ly./lx;        %interaction along y
    Dx = D./ly; %0.17;              %DM strength along x
    Dy = D./gamma.*(1./lx);         %DM strength along y
    K =  K./(gamma).*1./(lx.*ly);      %K > 0 is an easy-axis anisotropy
    
    %Computing the analytical texture
    [nx,ny,nz] = RelativisticProfile_tweaked(x,y,X,Y,beta,lx,ly,R,w,A);

    %Computing energy difference
    deltaE = abs(energysquare3(Sx,Sy,Sz,Jx,Jy,Dx,Dy,K) - energysquare3(nx,ny,nz,Jx,Jy,Dx,Dy,K));

end