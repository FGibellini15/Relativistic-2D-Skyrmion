function [R,w,A,Delta] = MinimizeDeltaE(Sx,Sy,Sz,J,K,beta,rho)
    %By minimizing the energy difference, the analytical profile
    %is fitted to the numerical skyrmion, denoted by Sx,Sy,Sz

    %Contraction coefficients and guesses
    [lx,ly] = ContractionCompute_raw(J,K,rho,beta);
    
    R_guess = 10;
    w_guess = 4.4;
    %lx_guess = 1;
    %ly_guess = 1;
    X_guess = [R_guess,w_guess,0];
   
    %Fitting
    [Y,Fmin] = fminsearch(@(X) energy_diff(Sx,Sy,Sz,J,K,beta,rho,lx,ly,X(1),X(2),X(3)), X_guess);
    
    % Extract the optimized parameters from the fitting result
    R = Y(1);
    w = Y(2);
    A = Y(3);
    
    Delta = Fmin;

end