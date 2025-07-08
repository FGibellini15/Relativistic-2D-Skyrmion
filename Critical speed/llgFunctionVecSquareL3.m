%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to calculate the LLG equation, the RHS of the differential eq:
% dS/dt = F(S,t)
%       = C(1)*SxHeff + C(2)*SxSxHeff
% where 
%       S     = spin vector
%       Heff  = efective field
%
% Heff is in turn a function of S, D and Hfield, where
%       J1 is assumed to be 1
%       D is the DM coupling constant
%       H is the applied magnetic field, vector of dimension 3
%
% This is the vector version, to be compatible with built-in ODE solvers
% This version is compatible with the square lattice configuration
% 
% A. P. - Groningen - Feb 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUT PARAMETERS:
% t, time (not used, but needed for ODE solver)
% Sx,Sy,Sz: matrices containing spin components at lattice sites (all in SS, column vector)
%+ SS details: SS is a column vector, made up of 3 matrices extended in columns
% C: vector of LLG coefficients
% D: the DM coupling constant
% H: H field applied [Hx Hy Hz]
%
function sDot = llgFunctionVecSquareL3(t,SS,C,Jx,Jy,Dx,Dy,K,Nx,Ny)
    N2 = Nx*Ny; % number of lattice points in one layer
    % umpack spin matrices to compute fields
    Sx = reshape(SS(1:N2),Ny,Nx);
    Sy = reshape(SS(N2+1:2*N2),Ny,Nx);
    Sz = reshape(SS(2*N2+1:3*N2),Ny,Nx);
    % compute the effective magnetic field for the bilayer
    [HeffX,HeffY,HeffZ] = EffFieldsquare3(Sx,Sy,Sz,Jx,Jy,Dx,Dy,K);
    % components of S x Heff cross product
    sXhX = Sy.*HeffZ-Sz.*HeffY;	% x component 
	sXhY = Sz.*HeffX-Sx.*HeffZ;	% y component 
	sXhZ = Sx.*HeffY-Sy.*HeffX;	% z component 
    % components of S x (S x Heff) cross product
    sXsXhX = Sy.*sXhZ - Sz.*sXhY;   % x component 
    sXsXhY = Sz.*sXhX - Sx.*sXhZ;   % y component 
    sXsXhZ = Sx.*sXhY - Sy.*sXhX;   % z component 
    % components of LLG eq.
    sDotX = C(1)*sXhX + C(2)*sXsXhX;
    sDotY = C(1)*sXhY + C(2)*sXsXhY;
    sDotZ = C(1)*sXhZ + C(2)*sXsXhZ;
    % pack them in a column vector
    sDot = [sDotX(1:N2) sDotY(1:N2) sDotZ(1:N2)]';
    %fprintf('CHECK: t = %f\n',t)
end
%EOF