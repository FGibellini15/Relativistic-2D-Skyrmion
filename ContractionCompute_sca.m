function [l_x,l_y] = ContractionCompute_sca(J,K,rho,beta)
    %This function computes the additional contraction
    %using also the scaling argument

    %Defines DM interaction coefficient
    D = rho.*2.*sqrt(K.*J)./pi;
    
    %Obtains static values for thickness and radius
    w = pi.*D./(2.*K);
    R = w./sqrt(1-rho.^2);
    
    %Obtains relativistic coefficient
    gamma = 1./(sqrt(1-beta.^2));
    
    %Evaluates integrals (R>>w)
    M = 4.*pi.*(R./w + w./R); %int of (grad n).^2
    N = 2.*(pi.^2).*R; %int of -(n . (curl of n))
    U = 4.*pi.*w.*R; %int of (1-nz.^2)
        
    %a,b,c as defined for the quadratic equation
    a = J .* M./2 .* (gamma.^2 - 1);
    b = (D.*N - 2.*J.*M).*gamma;
    c = 2.*J.*M - K.*U;

    %Obtaining the contraction coefficients
    l_x = (-b - sqrt(b.^2 - 4.*a.*c))./(2.*a);
    l_y = 2-gamma.*l_x;
    
end