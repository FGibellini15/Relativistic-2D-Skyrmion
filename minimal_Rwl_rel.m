function [min_Rw,ratio_Rw,min_lx,fmins] = minimal_Rwl_rel(beta_min,beta_max,N_beta,rho_min,...
    rho_max,N_rho,J,K,Rw_guess,lx_guess) 
    %This function computes R,w and lx such to minimize the approximate
    %(R>>w) formula for energy of the analytical skyrmion
    
    %Ranges of beta and rho
    beta_values = linspace(beta_min,beta_max,N_beta);
    rho_values = linspace(rho_min,rho_max,N_beta);

    %Allocating space
    min_Rw = zeros(2, N_beta, N_rho);
    min_lx = zeros(N_beta, N_rho);
    fmins = zeros(N_beta, N_rho);

    %Finding minimizing parameters and ground state energy
    for i = 1:N_beta
        for j = 1:N_rho
            [y,Fmin] = fminsearch(@(X) rel_energy1(X(3), beta_values(i), J,...
                K, rho_values(j), X(1), X(2)), [Rw_guess lx_guess]); %Minimizes H-Pv
            
            Fmin = Fmin + momentum(y(3),J,beta_values(i),y(1),y(2)); %Adding momentum to get total energy
            
            %Extracting values
            min_Rw(:, i,j) = y(1:2);
            min_lx(i,j) = y(3); 
            fmins(i,j) = Fmin;
        end
    end
   
    %Computing ratio (not used later)
    ratio_Rw = zeros(N_beta, N_rho);
    for i = 1:N_beta
        for j = 1:N_rho
            y = min_Rw(1,i,j)/min_Rw(2,i,j);
            ratio_Rw(i,j) = y;
        end
    end

end