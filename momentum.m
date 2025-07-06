function P = momentum(l_x,J,beta,R,w)
    %This function computes momentum of the skyrmion according
    %to the approximate solutions for the integrals

    gamma = 1/sqrt(1-beta^2);
    l_y = 2 - gamma*l_x;
    M = 4*pi*(R/w + w/R); %int of (grad n)^2
    P = M/2*J/2*(beta^2)*(gamma)*l_x/l_y;
end