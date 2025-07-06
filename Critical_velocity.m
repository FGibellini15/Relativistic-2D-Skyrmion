beta_min = 0.001; %skyrmion velocity
beta_max = 1;
N_beta = 10;
rho_min = 0.1;  %rho = pi*D/(2*sqrt(JK))
rho_max = 1;
N_rho = 10;
Rw_guess = [22 4.4]; %Guesses for minimizing parameters, taken from static case
lx_guess = 1;

beta_values = linspace(beta_min,beta_max,N_beta);
rho_values = linspace(rho_min,rho_max,N_rho);
[rrho,bb] = meshgrid(rho,beta);

J = 1;
K = 0.5;

[lx,ly] = ContractionCompute_sca(J,K,rrho,bb)

crit_rho = zeros(1,N_beta);


for i = 1:N_beta 
    for j = i:N_rho
        rho_pr = rho(j)/sqrt(lx(i,j)^2 + ly(i,j)^2);
        if rho_pr > 0.999
            crit_rho(i) = rho(j)
            i = N_rho;
        end
    end
end