function y = rel_energy1(l_x, beta,J,K,rho,R,w)
   
    %Computes the value of H-Pv given the contraction lx (= lambda_x in the
    % paper), the speed of the skyrmion beta, the interaction coefficients
    %J,K, rho (which allows to obtain D), and the skyrmion size parameters
    %R and w

    %Basic coefficients
    D = rho*2*sqrt(K*J)/pi;     %DM interactions
    gamma = 1/(sqrt(1-beta^2)); %Lorentz

    %integrals
    M = 4*pi*(R/w + w/R) ; %int of (grad n)^2
    N = 2*(pi^2)*R ; %int of -(n . (curl of n))
    U = 4*pi*w*R ; %int of (1-nz^2)
        
    %Deriving contraction in y
    l_y = 2 - gamma*l_x;
    %l_y = 2*K*V/(D*N) - gamma*l_x; %In case we don't want to employ
    %scaling
    
    %Coefficients modifying the interactions
    a_r = 1/(2*gamma) * (l_x/l_y + l_y/l_x);
    b_r = -1/(2*gamma) * (1/l_x + gamma/l_y);
    c_r = 1/(gamma) * (1/(l_x*l_y));

    y = a_r*J/2*M + b_r*D*N + c_r*K/2*U;

end

