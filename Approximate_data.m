close all
%This file obtains values of R,w,lx,ly to minimize the value of H-Pv
%under the approximation R >> w
%It should be run before Plotter.m

%Interaction coefficients
J = 1; %Heisenberg
K = 0.05; %Anisotropy

%Varied Parameters (range)
beta_min = 0.01; %skyrmion velocity
beta_max = 0.26;
N_beta = 101;
rho_min = 0.9;  %rho = pi*D/(2*sqrt(JK))
rho_max = 0.98;
N_rho = 101;
Rw_guess = [22 4.4]; %Guesses for minimizing parameters, taken from static case
lx_guess = 1;

%Finding minimizing parameters
[min_Rw,ratio_Rw,min_lx,fmins]=minimal_Rwl_rel(beta_min,beta_max,N_beta,rho_min,...
    rho_max,N_rho,J,K,Rw_guess,lx_guess);

%PLOTTING

%Defining spaces
beta_values = linspace(beta_min,beta_max,N_beta);
rho_values = linspace(rho_min,rho_max,N_rho);
[rrho,bb] = meshgrid(rho_values,beta_values);

R_array = reshape(min_Rw(1,:,:),[N_beta,N_rho]);
w_array = reshape(min_Rw(2,:,:),[N_beta,N_rho]);
lx_array = min_lx;

%Deriving ly
ly_array = zeros(size(lx_array)); %Obtaining l_y 
for i = 1:N_beta
    for j = 1:N_rho
        gamma = 1/sqrt(1-(bb(i,j))^2);
        ly_array(i,j) = 2 - gamma*lx_array(i,j);
    end
end

%l_x and l_y derived by analysis, for comparison
[lx_array_sca,ly_array_sca] = ContractionCompute_sca(J,K,rrho,bb);
[lx_array_raw,ly_array_raw] = ContractionCompute_raw(J,K,rrho,bb);

%Renaming of values for further purpose (Plotter.m)
E_ana = fmins;
R_ana = R_array;
w_ana = w_array;

rho_ana = rho_values;
beta_ana = beta_values;
rrho_ana = rrho;
bb_ana = bb;
N_r_ana = N_rho;
N_b_ana = N_beta;

%Useful values of beta and rho
rho1_a = rho_values(1);
rho2_a = rho_values(ceil(N_r_ana/2));
rho3_a = rho_values(N_r_ana);

beta1_a = beta_values(1);
beta2_a = beta_values(ceil(N_b_ana/2));
beta3_a = beta_values(N_b_ana);

rho1_a = round(rho1_a,2);
rho2_a = round(rho2_a,2);
rho3_a = round(rho3_a,2);

beta1_a = round(beta1_a,2);
beta2_a = round(beta2_a,2);
beta3_a = round(beta3_a,2);



% The actual plotting can be performed through the Plotter.m file
%In alternative, it can also be obtained by uncommenting the area below

% Plotting R
% Contours
% Fig_blR = figure; 
% contourf(rrho,bb, R_array, 100, LineStyle = 'None'); 
% c = colorbar  
% c.Label.String = 'asa';
% c.Label.Position(2) = c.Position(2) - 0.1; % Move label down
% xlabel('\rho',FontSize=20);
% ylabel('\beta',FontSize=20);
% title('Contour Plot of R vs \beta and \rho',FontSize=20);
% hold off
% 
% R in function of rho
% Fig_rR = figure;
% hold on
% plot(rho_values, R_array(1,:),DisplayName = ['\beta = ' num2str(beta1)]); 
% plot(rho_values, R_array(N_beta/2,:),DisplayName = ['\beta = ' num2str(beta2)]); %Need even N_beta
% plot(rho_values, R_array(N_beta,:),DisplayName = ['\beta = ' num2str(beta3)]); 
% 
% 
% xlabel('\rho',FontSize=20);
% ylabel('R',FontSize=20);
% title('Plot of R vs \rho',FontSize=20);
% legend('Location', 'northwest',FontSize=15);
% hold off
% 
% R in function of beta
% Fig_bR = figure;
% hold on
% plot(beta_values, R_array(:,1),DisplayName = ['\rho = ' num2str(rho1)]); 
% plot(beta_values, R_array(:,N_rho/2), DisplayName = ['\rho = ' num2str(rho2)]); 
% plot(beta_values, R_array(:,N_rho), DisplayName = ['\rho = ' num2str(rho3)]); 
% 
% axis tight;
% xlabel('\beta',FontSize=20);
% ylabel('R',FontSize=20);
% title('Plot of R vs \beta',FontSize=20);
% legend('Location', 'northwest',FontSize=15);
% hold off
% 
% Plotting w
% Contour
% Fig_blw = figure;
% contourf(rrho,bb, w_array, 100,LineStyle = 'None'); colorbar
% xlabel('\rho',FontSize=20);
% ylabel('\beta',FontSize=20);
% title('Contour Plot of w vs \beta and \rho',FontSize=20);
% 
% w in function of rho
% Fig_rw = figure;
% hold on
% plot(rho_values, w_array(1,:),DisplayName = ['\beta = ' num2str(beta1)]); 
% plot(rho_values, w_array(N_beta/2,:),DisplayName = ['\beta = ' num2str(beta2)]); %Need even N_beta
% plot(rho_values, w_array(N_beta,:),DisplayName = ['\beta = ' num2str(beta3)]); 
% 
% xlabel('\rho',FontSize=20);
% ylabel('w',FontSize=20);
% title('Plot of w vs \rho',FontSize=20);
% legend('Location', 'northwest',FontSize=15);
% hold off
% 
% w in function of beta
% Fig_bw = figure;
% hold on
% plot(beta_values, w_array(:,1),DisplayName = ['\rho = ' num2str(rho1)]); 
% plot(beta_values, w_array(:,N_rho/2),DisplayName = ['\rho = ' num2str(rho2)]); 
% plot(beta_values, w_array(:,N_rho),DisplayName = ['\rho = ' num2str(rho3)]); 
% 
% axis tight;
% xlabel('\beta',FontSize=20);
% ylabel('w',FontSize=20);
% title('Plot of w vs \beta',FontSize=20);
% legend('Location', 'northwest',FontSize=15);
% hold off
% 
% 
% Minimal energy vs rho and beta
% Fig_bre = figure;
% contourf(rrho,bb, fmins, 100, LineStyle = 'None'); colorbar
% xlabel('\rho',FontSize=20);
% ylabel('\beta',FontSize=20);
% title('Contour Plot of minimal H-Pv vs \beta and \rho',FontSize=20);
% hold off
% 
% Plotting l_x (and l_y)
% Contours
% Fig_brl = figure;
% contourf(rrho,bb, lx_array, 100, LineStyle = 'None'); colorbar
% xlabel('\rho',FontSize=20);
% ylabel('\beta',FontSize=20);
% title('Contour Plot of \lambda_x vs \beta and \rho',FontSize=20);
% hold off
% 
% l_x in function of rho
% Fig_rlx = figure;
% hold on
% plot(rho_values, lx_array(1,:),DisplayName = ['\beta = ' num2str(beta1)]); 
% plot(rho_values, lx_array(N_beta/2,:),DisplayName = ['\beta = ' num2str(beta2)]);
% plot(rho_values, lx_array(N_beta,:),DisplayName = ['\beta = ' num2str(beta3)]); 
% plot(rho_values, lx_array_raw(1,:),':',DisplayName = ['\beta = ' num2str(beta1)]);
% plot(rho_values, lx_array_raw(N_beta/2,:),':',DisplayName = ['\beta = ' num2str(beta2)]);
% plot(rho_values, lx_array_raw(N_beta,:),':',DisplayName = ['\beta = ' num2str(beta3)]); 
% 
% 
% xlabel('\rho',FontSize=20);
% ylabel('\lambda_x',FontSize=20);
% title('Plot of \lambda_x vs \rho',FontSize=20);
% legend('Location', 'southeast',FontSize=15);
% hold off
% 
% l_x in function of beta
% Fig_blx = figure;
% hold on
% plot(beta_values, lx_array(:,1),DisplayName = ['\rho = ' num2str(rho1)]); 
% plot(beta_values, lx_array(:,N_rho/2),DisplayName = ['\rho = ' num2str(rho2)]); 
% plot(beta_values, lx_array(:,N_rho),DisplayName = ['\rho = ' num2str(rho3)]); 
% plot(beta_values, lx_array_raw(:,1),':', DisplayName = ['\rho = ' num2str(rho1)]); 
% plot(beta_values, lx_array_raw(:,N_rho/2),':',DisplayName = ['\rho = ' num2str(rho2)]); 
% plot(beta_values, lx_array_raw(:,N_rho),':',DisplayName = ['\rho = ' num2str(rho3)]); 
% 
% axis tight;
% xlabel('\beta',FontSize=20);
% ylabel('\lambda_x',FontSize=20);
% title('Plot of \lambda_x vs \beta',FontSize=20);
% legend('Location', 'southwest',FontSize=15);
% hold off
% 
% l_y in function of rho
% Fig_rly = figure;
% hold on
% plot(rho_values, ly_array(1,:),DisplayName = ['\beta = ' num2str(beta1)]); 
% plot(rho_values, ly_array(N_beta/2,:),DisplayName = ['\beta = ' num2str(beta2)]);
% plot(rho_values, ly_array(N_beta,:),DisplayName = ['\beta = ' num2str(beta3)]); 
% plot(rho_values, ly_array_raw(1,:),':',DisplayName = ['\beta = ' num2str(beta1)]); 
% plot(rho_values, ly_array_raw(N_beta/2,:),':',DisplayName = ['\beta = ' num2str(beta2)]);
% plot(rho_values, ly_array_raw(N_beta,:),':',DisplayName = ['\beta = ' num2str(beta3)]); 
% 
% axis tight;
% xlabel('\rho',FontSize=20);
% ylabel('\lambda_y',FontSize=20);
% title('Plot of \lambda_y vs \rho',FontSize=20);
% legend('Location', 'southwest',FontSize=15);
% hold off
% 
% l_y in function of beta
% Fig_bly = figure;
% hold on
% plot(beta_values, ly_array(:,1),DisplayName = ['\rho = ' num2str(rho1)]); 
% plot(beta_values, ly_array(:,N_rho/2),DisplayName = ['\rho = ' num2str(rho2)]); 
% plot(beta_values, ly_array(:,N_rho),DisplayName = ['\rho = ' num2str(rho3)]); 
% plot(beta_values, ly_array_raw(:,1),':', DisplayName = ['\rho = ' num2str(rho1)]); 
% plot(beta_values, ly_array_raw(:,N_rho/2),':',DisplayName = ['\rho = ' num2str(rho2)]); 
% plot(beta_values, ly_array_raw(:,N_rho),':',DisplayName = ['\rho = ' num2str(rho3)]); 
% 
% axis tight;
% xlabel('\beta',FontSize=20);
% ylabel('\lambda_y',FontSize=20);
% title('Plot of \lambda_y vs \beta',FontSize=20);
% legend('Location', 'southwest',FontSize=15);
% hold off
% 
% %R in function of l_x
% Fig_lRb = figure;
% hold on
% plot(lx_array(:,1), R_array(:,1),DisplayName = DisplayName = ['\rho = ' num2str(rho1)]); 
% plot(lx_array(:,N_rho/2), R_array(:,N_rho/2),DisplayName = ['\rho = ' num2str(rho2)]); 
% plot(lx_array(:,N_rho), R_array(:,N_rho),DisplayName = ['\rho = ' num2str(rho3)]); 
% 
% xlabel('l_x');
% ylabel('R');
% title('Plot of R vs l_x at constant beta');
% legend('Location', 'northwest',FontSize=15);
% hold off
% 
% Fig_lRrho = figure;
% hold on
% plot(lx_array(1,:), R_array(1,:),DisplayName = ['\beta = ' num2str(beta1)]); 
% plot(lx_array(N_beta/2,:), R_array(N_beta/2,:),DisplayName = ['\beta = ' num2str(beta2)]); 
% plot(lx_array(N_beta,:), R_array(N_beta,:),DisplayName = ['\beta = ' num2str(beta3)]); 
% 
% xlabel('l_x');
% ylabel('R');
% title('Plot of R vs l_x at constant rho');
% legend('Location', 'northwest',FontSize=15);
% hold off
% 
% %w in function of l_x
% Fig_lwb = figure;
% hold on
% plot(lx_array(:,1), w_array(:,1),DisplayName = DisplayName = ['\rho = ' num2str(rho1)]); 
% plot(lx_array(:,N_rho/2), w_array(:,N_rho/2),DisplayName = ['\rho = ' num2str(rho2)]); 
% plot(lx_array(:,N_rho), w_array(:,N_rho),DisplayName = ['\beta = ' num2str(beta3)]); 
% 
% xlabel('l_x');
% ylabel('w');
% title('Plot of w vs l_x at constant beta');
% legend('Location', 'northwest',FontSize=15);
% hold off
% 
% Fig_lwrho = figure;
% hold on
% plot(lx_array(1,:), w_array(1,:),DisplayName = ['\beta = ' num2str(beta1)]); 
% plot(lx_array(N_beta/2,:), w_array(N_beta/2,:),DisplayName = ['\beta = ' num2str(beta2)]); 
% plot(lx_array(N_beta,:), w_array(N_beta,:),DisplayName = ['\beta = ' num2str(beta3)]); 
% 
% xlabel('l_x');
% ylabel('w');
% title('Plot of w vs l_x at constant rho');
% legend('Location', 'northwest',FontSize=15);
% hold off
