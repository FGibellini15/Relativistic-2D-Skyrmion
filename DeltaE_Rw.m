function [E0,E0_model,deltaE,R,w] = DeltaE_Rw(Sx,Sy,Sz,J,K,beta,rho,time,time_step,aRel)
    
    %This function computes takes as input a numerical texture and relaxes
    %it, under the given interaction coefficients J,K, rho and beta
    %(skyrmion velocity). A skyrmion witht the analytical profile is then
    %fitted to the numerical relaxed one. It returns the radius and thickness
    % so obtained, together with the energy difference. 
    % The function computes and returns also the total energy of both configurations.
    
    %Generates lattice for the fitted skyrmion
    [Ny,Nx] = size(Sx); %size, must MAKE Nx and Ny ODD
    [x,y] = meshgrid(0:1:Nx-1,0:1:Ny-1);
    X = floor((Nx-1)/2); %Coordinates of the center
    Y = floor((Ny-1)/2);

    %Basic and Relativistic coefficients
    gamma = 1./sqrt(1-beta.^2);
    D = rho*2*sqrt(K*J)/pi;
    [lx,ly] = ContractionCompute_raw(J,K,rho,beta);
    
    %Interaction Coefficients
    Jx = J./(gamma).*lx./ly;         %interaction along x
    Jy = J./(gamma).*ly./lx;        %interaction along y
    Dx = D./ly; %0.17;   %DM strength along x
    Dy = D./gamma.*(1./lx);         %DM strength along y
    K =  K./(gamma).*1./(lx.*ly);      %K > 0 is an easy-axis anisotropy
   
    %Relaxing given skyrmion
    [Sx,Sy,Sz] = Relax(Sx,Sy,Sz,Jx,Jy,Dx,Dy,K,time,time_step,aRel);

    %Computing minimizing parameters
    [R,w,A,deltaE] = MinimizeDeltaE(Sx,Sy,Sz,J,K,beta,rho);
    [nx,ny,nz] = RelativisticProfile_tweaked(x,y,X,Y,beta,lx,ly,R,w,A);       
    
    %Computing Total Energy (adding momentum) 
    Jx = Jx*gamma^2*(1+beta^2);
    E0 = energysquare3(Sx,Sy,Sz,Jx,Jy,Dx,Dy,K);
    E0_model = energysquare3(nx,ny,nz,Jx,Jy,Dx,Dy,K);
    
end