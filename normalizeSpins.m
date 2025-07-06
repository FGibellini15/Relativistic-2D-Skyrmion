function [Nx,Ny,Nz] = normalizeSpins(Sx,Sy,Sz)
    % Normalizes a given texture
    SS = sqrt(Sx.^2+Sy.^2+Sz.^2); 
    Nx = Sx./SS;
    Ny = Sy./SS;
    Nz = Sz./SS;  
end
%EOF