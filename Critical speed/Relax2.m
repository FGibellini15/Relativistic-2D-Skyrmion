function [Sxc,Syc,Szc,E] = Relax2(Sx,Sy,Sz,Jx,Jy,Dx,Dy,K,time,time_step,aRel)

%Takes as argument the spin texture to be relaxed, the interaction coefficients 
%in the single directions and the time and damping parameters for the
%relaxation
[Ny,Nx] = size(Sx); %size, MAKE Nx and Ny ODD!!!
NSITES = Nx*Ny;     %number of sites
Kx = floor((Nx - 1)/2);
Ky = floor((Ny - 1)/2);
x = (0:Nx-1) - Kx;
y = (0:Ny-1) - Ky;
[X,Y] = meshgrid(x,y);
%%%%%%%%%%%%%%%%%%%%%%
%relaxation parameters
%%%%%%%%%%%%%%%%%%%%%%
tStart = 0;             % start time 
tEnd = time;            % finish time. EST ex. time per unit time: 0.00131314
NM = time_step;                % time interval between measurements
tSpan = tStart:NM:tEnd; % time vector for measurements

%LLG equation coefficients. NB: using signs from DOI: 10.1038/NMAT3862

Crel=zeros(1,2);        % initialize vector of constants
Crel(1)=-1/(1+aRel^2);  % coeff. for SxHeff
Crel(2)=aRel*Crel(1);   % coeff. for SxSxHeff
%%%%%%%%%%% 
%simulation
%%%%%%%%%%%
%Initial spin configuration
%normalize spins
%[Sx,Sy,Sz] = normalizeSpins(Sx,Sy,Sz);
%plot

pause(3)
% save initial spin configuration
Sx0 = Sx;
Sy0 = Sy;
Sz0 = Sz;
%%INITIALIZE VECTORS
Nt = length(tSpan);    %number of time points measured
Et = zeros(1,Nt);      %energy as a function of time
%% INTEGRATION
%Initial conditions for LLG
%pack the spin matrices into one big column vector (to use with built-in ODE solver)
SS0 = [Sx(1:NSITES) Sy(1:NSITES) Sz(1:NSITES)]';
%
%Initial energy
EFM = -(Jx+Jy); %energy of FM state per site, S perpendicular to z
E0 = energysquare3(Sx,Sy,Sz,Jx,Jy,Dx,Dy,K) - NSITES*EFM;
%solve ode
%---------
tic
    options = odeset('RelTol',1e-6,'AbsTol',1e-10,'InitialStep',1e-3,'MaxStep',1.,'stats','on');
    [T,ST] = ode45(@(t,SS) llgFunctionVecSquareL3(t,SS,Crel,Jx,Jy,Dx,Dy,K,Nx,Ny),tSpan,SS0,options);
    %T: Column vector of time points.
    %ST: Solution array. Each row in ST corresponds to the solution at a time returned in the corresponding row of T.
toc
% number of solution time points
%[nT,~] = size(ST);
%% ANALYSIS 
% analyze the S matrices and compute interesting quantities
for mit=1:Nt
 %spin configuration
 SSc = ST(mit,:); %vector form
 % unpack spin matrices of final configuration

 NX = Nx;
 NY = Ny;
Sxc = reshape(SSc(1:NSITES),NY,NX);
Syc = reshape(SSc(NSITES+1:2*NSITES),NY,NX);
Szc = reshape(SSc(2*NSITES+1:3*NSITES),NY,NX);
% normalize the solution
 [Sxc,Syc,Szc] = normalizeSpins(Sxc,Syc,Szc);
 %energy
 Etmp = energysquare3(Sxc,Syc,Szc,Jx,Jy,Dx,Dy,K) - NSITES*EFM;
 Et(mit) = Etmp;
end
% final spin configuration
Sx = Sxc;
Sy = Syc;
Sz = Szc;

E = Et(Nt);

%----------------------
%energy plot
%----------------------
% fprintf('Initial Energy:          E = %f\n',E0);
% fprintf('Energy after relaxation: E = %f\n',Et(Nt));
% Einitial = E0;
% Efinal = Et(Nt);
% he = plot(T,Et,'ro');
% set(gca,'Fontsize',16)
% xlabel(' Time ')
% ylabel(' Energy ')
% pause(3)
%saving relaxed spin configuration
% Conical Spiral
%save meronsK0_01dJ0_25J20_3.mat Sx Sy Sz deltaJ J2 J3 K HField Efinal 
%
%figure
end