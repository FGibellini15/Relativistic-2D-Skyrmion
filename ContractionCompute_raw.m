function [l_x,l_y] = ContractionCompute_raw(J,K,rho,beta)
    %This function computes the additional contraction without
    %resorting to the scaling argument

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
        
    %Defining useful coefficient
    alpha = (2.*K.*U)./(D.*N);
    
    %a,b,c as defined for the quadratic equation
    a = J .* M./2 .* (gamma.^2 - 1);
    b = (D.*N - alpha.*J.*M).*gamma;
    c = J.*M.*(alpha.^2)./2 - K.*U;
    
    %Obtaining the contraction coefficients
    l_x = (-b - sqrt(b.^2 - 4.*a.*c))./(2.*a);
    l_y = alpha - gamma.*l_x; 

end