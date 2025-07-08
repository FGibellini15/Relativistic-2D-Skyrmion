close all
%This file computes the critical speed of the numerical skyrmion

%Defining ranges
beta_min = 0.6;
beta_max = 0.9;
N_beta = 3;
beta_values = linspace(beta_min, beta_max, N_beta);
gamma_values = 1./sqrt(1-beta_values.^2);
rho_min = 0.1;
rho_max = 0.99;
N_rho = 10;
rho_values = linspace(rho_min, rho_max, N_rho);

%Allocating space
v_crit = ones(N_rho,1);
rho_crit = ones(N_beta,1);

E = zeros(N_rho,N_beta);

for i = 1:N_beta
    for j=7
        load skyrmion5_71_min1.mat;
        Jx = 1.0;         %interaction along x
        Jy = 1.0;         %interaction along y
        J = 1.0;
        D = rho_values(j)*2*sqrt(K*J)/pi;     %DM interactions (from rho)
        Dy = D;         %DM strength along y
        K =  0.055;      %K > 0 is an easy-axis anisotropy

        Dx = D*(gamma_values(i)); %DM strength along x, adjusted for gamma

        [Sx,Sy,Sz,E(j,i)] = Relax2(Sx,Sy,Sz,Jx,Jy,Dx,Dy,K,800,200,1);
        
        tit = ['\rho = ' num2str(round(rho_values(j),2)) ' \beta = ' num2str(round(beta_values(i),2))];
        Texture_plot_contr(71,71,Sx,Sy,Sz,gamma_values(i),tit)
        
        % if i>1 && abs(E(j,i) - E(j,i-1))>0.3*E(j,i-1)
        % 
        %    v_crit(j) = beta_values(i-1) 
        % 
        %    break;
        % end

        %Texture_plot(71,71,Sx,Sy,Sz);
        %E = energysquare3(Sx,Sy,Sz,Jx,Jy,Dx,Dy,K)
    end
end




% for i = N_beta
%     for j=1:N_rho
%         load skyrmion5_71_min1.mat;
%         Jx = 1.0;         %interaction along x
%         Jy = 1.0;         %interaction along y
%         D = rho_values(j)*2*sqrt(K*J)/pi;     %DM interactions (from rho)
%         Dy = D;         %DM strength along y
%         K =  0.055;      %K > 0 is an easy-axis anisotropy
% 
%         Dx = D*(gamma_values(i)); %DM strength along x, adjusted for gamma
% 
%         [Sx,Sy,Sz,E(j,i)] = Relax2(Sx,Sy,Sz,Jx,Jy,Dx,Dy,K,800,200,1);
%         % 
%         if E(j,i)<0
% 
%            rho_crit(i) = rho_values(j) 
% 
%            break;
%         end
% 
%         %Texture_plot(71,71,Sx,Sy,Sz);
%         %E = energysquare3(Sx,Sy,Sz,Jx,Jy,Dx,Dy,K)
%     end
% end
% 
% 
% 



% Fig = figure;
% scatter(rho_values,E(:,7))
% % Fig = figure;
% % plot(beta_values,E(3,:))
% % Fig = figure;
% % plot(beta_values,E(5,:))
% Fig = figure;
% plot(beta_values,rho_crit)
