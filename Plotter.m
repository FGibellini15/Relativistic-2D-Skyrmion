%This file plots all the figures presented in the paper.
%It should be run after Numerical_and_fitted_data.m and Approximate_data.m

close all

%%PLOTTING
color1 = "#30a8c6";
color2 = "#b53f26";
color3 = "#097969";

%% Plotting R
%Contours
Fig_1 = figure; 

subplot(1,2,1)

contourf(rrho_ana,bb_ana, R_ana, 100, LineStyle = 'None'); 
c = colorbar;  
c.FontSize = 12;
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
xlabel('\rho',FontSize=20);
ylabel('\beta',FontSize=20);
title('(a) Analytical R ',FontSize=20);
hold off

subplot(1,2,2)
contourf(rrho_fit,bb_fit, R_fit, 100, LineStyle = 'None'); 
c = colorbar;  
c.FontSize = 12;
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
xlabel('\rho',FontSize=20);
ylabel('\beta',FontSize=20);
title('(b) Numerical R',FontSize=20);
hold off

%R in function of rho
Fig2 = figure;
%sgtitle('Plot of R vs \rho',FontSize=20);

subplot(1,3,1)
hold on
plot(rho_fit, R_fit(1,:),DisplayName = 'Fitted', Color=color1); 
plot(rho_ana, R_ana(1,:),DisplayName = 'Analytical',LineStyle = '--', LineWidth=2, Color=color1); 

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
title('(a)', FontSize = 15)

xlabel('\rho',FontSize=20);
ylabel('R',FontSize=20);
leg = legend('Location', 'northwest',FontSize=15);
title(leg,['\beta = ' num2str(beta1_fit)])

hold off

subplot(1,3,2)
hold on
plot(rho_fit, R_fit(ceil(N_b_fit/2),:),DisplayName = 'Fitted', Color=color2); %Need even N_b_fit
plot(rho_ana, R_ana(ceil(N_b_ana/2),:),DisplayName = 'Analytical',LineStyle = '--', LineWidth=2, Color=color2); %Need even N_b_fit

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
title('(b)', FontSize = 15)

xlabel('\rho',FontSize=20);
%ylabel('R',FontSize=20);
leg = legend('Location', 'northwest',FontSize=15);
title(leg,['\beta = ' num2str(beta2_fit)])

hold off

subplot(1,3,3)
hold on
plot(rho_fit, R_fit(N_b_fit,:),DisplayName = 'Fitted', Color=color3); 
plot(rho_ana, R_ana(N_b_ana,:),DisplayName = 'Analytical',LineStyle = '--', LineWidth=2, Color=color3); 

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
title('(c)', FontSize = 15)

leg = legend('Location', 'northwest',FontSize=15);
title(leg,['\beta = ' num2str(beta3_fit)])
xlabel('\rho',FontSize=20);
%ylabel('R',FontSize=20);
hold off


%R in function of beta
Fig3 = figure;
%sgtitle('Plot of R vs \beta',FontSize=20);

subplot(1,3,1);
hold on
plot(beta_fit, R_fit(:,1),DisplayName = 'Fitted', Color=color1); 
plot(beta_ana, R_ana(:,1),DisplayName = 'Analytical',LineStyle = '--', LineWidth=2, Color=color1); 

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
title('(a)', FontSize = 15)

xlabel('\beta',FontSize=20);
ylabel('R',FontSize=20);
l = legend('Location', 'northwest',FontSize=15);
title(l,['\rho = ' num2str(rho1_fit)])

bottom = R_ana(1,1)-0.5;
axis tight;
ylim([bottom inf])

hold off


subplot(1,3,2);
hold on
plot(beta_fit, R_fit(:,ceil(N_r_fit/2)),DisplayName = 'Fitted', Color=color2); 
plot(beta_ana, R_ana(:,ceil(N_r_ana/2)), DisplayName = 'Analytical', LineStyle = '--', LineWidth=2,Color=color2); 

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
title('(b)', FontSize = 15)

xlabel('\beta',FontSize=20);
%ylabel('R',FontSize=20);
l = legend('Location', 'northwest',FontSize=15);
title(l,['\rho = ' num2str(rho2_fit)]);

bottom = R_ana(1,ceil(N_r_ana/2))-0.5;
axis tight;
ylim([bottom inf])

hold off


subplot(1,3,3);
hold on
plot(beta_fit, R_fit(:,N_r_fit), DisplayName = 'Fitted', Color=color3); 
plot(beta_ana, R_ana(:,N_r_ana), DisplayName = 'Analytical',LineStyle = '--', LineWidth=2,Color=color3); 

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
title('(c)', FontSize = 15)
bottom = R_fit(1,N_r_fit)-1;
axis tight;
ylim([bottom 50])

xlabel('\beta',FontSize=20);
%ylabel('R',FontSize=20);
legend('Location', 'northwest',FontSize=15);
l = legend('Location', 'northwest',FontSize=15);
title(l,['\rho = ' num2str(rho3_fit)])


hold off



%% Plotting w
%Contour
Fig_blwa = figure;

subplot(1,2,1)
contourf(rrho_ana,bb_ana, w_ana, 100,LineStyle = 'None'); 
c = colorbar;

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
c.FontSize = 12;

xlabel('\rho',FontSize=20);
ylabel('\beta',FontSize=20);
title('(a) Analytical w',FontSize=20);

%Fig_blwf = figure;
subplot(1,2,2)
contourf(rrho_fit,bb_fit, w_fit, 100,LineStyle = 'None');
c = colorbar;

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
c.FontSize = 12;

xlabel('\rho',FontSize=20);
ylabel('\beta',FontSize=20);
title('(b) Numerical w',FontSize=20);


%w in function of rho
Fig_rw = figure;
%sgtitle('Plot of w vs \rho',FontSize=20);

subplot(1,3,1);
hold on
plot(rho_fit, w_fit(1,:),DisplayName = 'Fitted', Color=color1); 
plot(rho_ana, w_ana(1,:),DisplayName = 'Analytical', LineStyle = '--', LineWidth=2, Color = color1); 

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
title('(a)', FontSize = 15)

xlabel('\rho',FontSize=20);
ylabel('w',FontSize=20);
l = legend('Location', 'northwest',FontSize=15);
title(l,['\beta = ' num2str(beta1_fit)])
hold off

subplot(1,3,2)
hold on
plot(rho_fit, w_fit(ceil(N_b_fit/2),:),DisplayName = 'Fitted', Color=color2); %Need even N_b_fit
plot(rho_ana, w_ana(ceil(N_b_ana/2),:),DisplayName = 'Analytical', LineStyle = '--', LineWidth=2, Color = color2); %Need even N_b_fit

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
title('(b)', FontSize = 15)

xlabel('\rho',FontSize=20);
%ylabel('w',FontSize=20);
legend('Location', 'northwest',FontSize=15);
l = legend('Location', 'northwest',FontSize=15);
title(l,['\beta = ' num2str(beta2_fit)])

hold off

subplot(1,3,3)
hold on
plot(rho_fit, w_fit(N_b_fit,:),DisplayName = 'Fitted', Color=color3); 
plot(rho_ana, w_ana(N_b_ana,:),DisplayName = 'Analytical', LineStyle = '--', LineWidth=2,Color = color3); 

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
title('(c)', FontSize = 15)

xlabel('\rho',FontSize=20);
%ylabel('w',FontSize=20);
legend('Location', 'northwest',FontSize=15);
l = legend('Location', 'northwest',FontSize=15);
title(l,['\beta = ' num2str(beta3_fit)])

hold off

%w in function of beta
Fig_bw = figure;
%sgtitle('Plot of w vs \beta',FontSize=20);

subplot(1,3,1)
hold on
plot(beta_fit, w_fit(:,1),DisplayName = 'Fitted',Color = color1); 
plot(beta_ana, w_ana(:,1),DisplayName = 'Analytical',LineStyle = '--', LineWidth=2,Color = color1); 

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
title('(a)', FontSize = 15)
axis tight;


xlabel('\beta',FontSize=20);
ylabel('w',FontSize=20);
l = legend('Location', 'northeast',FontSize=15);
title(l,['\rho = ' num2str(rho1_fit)]);
hold off

subplot(1,3,2)
hold on
plot(beta_fit, w_fit(:,ceil(N_r_fit/2)),DisplayName = 'Fitted',Color = color2); 
plot(beta_ana, w_ana(:,ceil(N_r_ana/2)),DisplayName = 'Analytical', LineStyle = '--', LineWidth=2,Color = color2); 

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
title('(b)', FontSize = 15)
top = w_ana(1,ceil(N_r_ana/2))+0.02;
axis tight;
ylim([-inf top])

xlabel('\beta',FontSize=20);
%ylabel('w',FontSize=20);
l = legend('Location', 'northwest',FontSize=15);
title(l,['\rho = ' num2str(rho2_fit)]);
hold off

subplot(1,3,3)
hold on
plot(beta_fit, w_fit(:,N_r_fit),DisplayName = 'Fitted',Color = color3); 
plot(beta_ana, w_ana(:,N_r_ana),DisplayName = 'Analytical', LineStyle = '--', LineWidth=2,Color = color3); 


ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
title('(c)', FontSize = 15)
top = w_ana(1,N_r_ana)+0.02;
axis tight;
ylim([-inf top])

xlabel('\beta',FontSize=20);
%ylabel('w',FontSize=20);
l = legend('Location', 'northwest',FontSize=15);
title(l,['\rho = ' num2str(rho3_fit)]);
hold off


%% Plotting Delta E
%Minimal deltaE vs rho and beta
Fig_brde = figure;
sgtitle('Minimal \Delta(H-Pv) vs \beta and \rho',FontSize=20);

subplot(1,3,1);
hold on
contourf(rrho_fit,bb_fit, deltaE_fit, 100, LineStyle = 'None'); 
c = colorbar;

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
c.FontSize = 12;
title('(a)',FontSize = 15);

xlabel('\rho',FontSize=20);
ylabel('\beta',FontSize=20);
hold off

%DeltaE vs beta
subplot(1,3,2);
hold on
plot(beta_fit, deltaE_fit(:,1),DisplayName = [num2str(rho1_fit)],Color = color1); 
plot(beta_fit, deltaE_fit(:,ceil(N_r_fit/2)),DisplayName = [num2str(rho2_fit)],Color = color2); 
plot(beta_fit, deltaE_fit(:,N_r_fit),DisplayName = [num2str(rho3_fit)],Color = color3); 

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
title('(b)', FontSize = 15)
bottom = deltaE_fit(1,1)-0.02;
axis tight;
ylim([bottom inf])

xlabel('\beta',FontSize=20);
ylabel('\Delta(H-Pv)',FontSize=20);
%title('Plot of \Delta E vs \beta',FontSize=20);
l = legend('Location', 'northwest',FontSize=15);
title(l,'\rho');
hold off

%DeltaE vs rho
subplot(1,3,3);
hold on
plot(rho_fit, deltaE_fit(1,:),DisplayName = [num2str(rho1_fit)],Color = color1); 
plot(rho_fit, deltaE_fit(ceil(N_b_fit/2),:),DisplayName = [num2str(rho2_fit)],Color = color2); 
plot(rho_fit, deltaE_fit(N_b_fit,:),DisplayName = [num2str(rho3_fit)],Color = color3); 

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
title('(c)', FontSize = 15)
bottom = deltaE_fit(1,1)-0.02;
axis tight;
ylim([bottom inf])

xlabel('\rho',FontSize=20);
ylabel('\Delta(H-Pv)',FontSize=20);
%title('Plot of \Delta E vs \rho',FontSize=20);
l = legend('Location', 'northwest',FontSize=15);
title(l,'\rho');
hold off

%% Plotting Energy

gamma_ana = 1./(sqrt(1-beta_ana.^2));
gamma_fit = 1./(sqrt(1-beta_fit.^2));

%Analytical E vs gamma
Fig_bE = figure;
%sgtitle('Variation of H against \beta',FontSize=20);

subplot(1,3,1)
hold on
 
plot(gamma_fit, E_fit(:,1),DisplayName = 'Fitted', Color = color1); 
plot(gamma_ana, E_ana(:,1),DisplayName = 'Analytical' ,Color = color2);
plot(gamma_fit, gamma_fit*E_fit(1,1),DisplayName = 'Lor Fit',LineStyle = '--',LineWidth=2, Color = color1); 
plot(gamma_ana, gamma_ana*E_ana(1,1),DisplayName = 'Lor Ana',LineStyle = '--',LineWidth=2,Color = color2); 

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
title('(a)', FontSize = 15)
bottom = E_ana(1,1)-0.1;
axis tight;
ylim([bottom inf])

xlabel('\gamma',FontSize=20);
ylabel('E',FontSize=20);
l = legend('Location', 'northwest',FontSize=15);
title(l,['\rho = ' num2str(rho1_fit)]);
hold off

subplot(1,3,2)
hold on
plot(gamma_fit, E_fit(:,ceil(N_r_fit/2)),DisplayName = 'Fitted', Color = color1); 
plot(gamma_ana, E_ana(:,ceil(N_r_ana/2)),DisplayName = 'Analytical', Color = color2); 
plot(gamma_fit, gamma_fit*E_fit(1,ceil(N_r_fit/2)),DisplayName = 'Lor Fit', LineStyle='--', LineWidth=2, Color = color1); 
plot(gamma_ana, gamma_ana*E_ana(2,ceil(N_r_ana/2)),DisplayName = 'Lor Ana',LineStyle='--' , LineWidth=2, Color = color2); 

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
title('(b)', FontSize = 15)
bottom = E_ana(1,ceil(N_r_ana/2))-0.1;
axis tight;
ylim([bottom inf])

xlabel('\gamma',FontSize=20);
%ylabel('E',FontSize=20);
l = legend('Location', 'northwest',FontSize=15);
title(l,['\rho = ' num2str(rho2_fit)]);
hold off

subplot(1,3,3)
hold on
plot(gamma_fit, E_fit(:,N_r_fit),DisplayName = 'Fitted', Color = color1); 
plot(gamma_ana, E_ana(:,N_r_ana),DisplayName = 'Analytical', Color = color2); 
plot(gamma_fit, gamma_fit*E_fit(1,N_r_fit),DisplayName = 'Lor Fit', LineStyle='--', LineWidth=2, Color = color1); 
plot(gamma_ana, gamma_ana*E_ana(2,N_r_ana),DisplayName = 'Lor Ana',LineStyle='--' , LineWidth=2, Color = color2); 

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
title('(c)', FontSize = 15)
bottom = E_ana(1,N_r_ana)-0.1;
axis tight;
ylim([bottom inf])

xlabel('\gamma',FontSize=20);
%ylabel('E',FontSize=20);
l = legend('Location', 'northwest',FontSize=15);
title(l,['\rho = ' num2str(rho3_fit)]);
hold off


%E vs beta
Fig_bE = figure;
%sgtitle('Variation of H against \beta',FontSize=20);

subplot(1,3,1)
hold on
plot(beta_fit, E_fit(:,1),DisplayName = 'Fitted',Color = color1); 
plot(beta_fit, E_num(:,1),DisplayName = 'Numerical', LineStyle=":",LineWidth=2,Color = color1); 
plot(beta_ana, E_ana(:,1),DisplayName = 'Analytical', LineStyle = '--', LineWidth=2,Color = color1); 

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
title('(a)', FontSize = 15)
bottom = E_ana(1,1)-0.1;
axis tight;
ylim([bottom inf])

xlabel('\beta',FontSize=20);
ylabel('E',FontSize=20);
l = legend('Location', 'northwest',FontSize=15);
title(l,['\rho = ' num2str(rho1_fit)]);
hold off

subplot(1,3,2)
hold on
plot(beta_fit, E_fit(:,ceil(N_r_fit/2)),DisplayName = 'Fitted',Color = color2); 
plot(beta_fit, E_num(:,ceil(N_r_fit/2)),DisplayName = 'Numerical', LineStyle=":",LineWidth=2,Color = color2); 
plot(beta_ana, E_ana(:,ceil(N_r_ana/2)),DisplayName = 'Analytical', LineStyle='--', LineWidth=2, Color = color2); 

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
title('(b)', FontSize = 15)
bottom = E_ana(1,ceil(N_r_ana/2))-0.1;
axis tight;
ylim([bottom inf])

xlabel('\beta',FontSize=20);
%ylabel('E',FontSize=20);
l = legend('Location', 'northwest',FontSize=15);
title(l,['\rho = ' num2str(rho2_fit)]);
hold off

subplot(1,3,3)
hold on
plot(beta_fit, E_fit(:,N_r_fit),DisplayName = 'Fitted',Color = color3); 
plot(beta_fit, E_num(:,N_r_fit),DisplayName = 'Numerical', LineStyle=":",LineWidth=2,Color = color3); 
plot(beta_ana, E_ana(:,N_r_ana),DisplayName = 'Analytical', LineStyle='--', LineWidth=2, Color = color3); 

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
title('(c)', FontSize = 15)
bottom = E_ana(1,N_r_ana)-0.1;
axis tight;
ylim([bottom inf])

xlabel('\beta',FontSize=20);
%ylabel('E',FontSize=20);
l = legend('Location', 'northwest',FontSize=15);
title(l,['\rho = ' num2str(rho3_fit)]);
hold off

%E vs rho
Fig_rE = figure;
%sgtitle('Variation of H against \rho',FontSize=20);

subplot(1,3,1)
hold on
plot(rho_fit, E_fit(1,:),DisplayName = 'Fitted',Color = color1); 
plot(rho_fit, E_num(1,:),DisplayName = 'Numerical',LineStyle=":",LineWidth=2,Color = color1); 
plot(rho_ana, E_ana(1,:),DisplayName = 'Analytical',LineStyle = '--', LineWidth=2,Color = color1); 

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
title('(a)', FontSize = 15)
axis tight;

xlabel('\rho',FontSize=20);
ylabel('E',FontSize=20);
l = legend('Location', 'northeast',FontSize=15);
title(l,['\beta = ' num2str(beta1_fit)]);
hold off

subplot(1,3,2)
hold on
plot(rho_fit, E_fit(ceil(N_b_fit/2),:),DisplayName = 'Fitted',Color = color2); 
plot(rho_fit, E_num(ceil(N_b_fit/2),:),DisplayName = 'Numerical',LineStyle=":",LineWidth=2,Color = color2); 
plot(rho_ana, E_ana(ceil(N_b_ana/2),:),DisplayName = 'Analytical',LineStyle = '--', LineWidth=2,Color = color2); 

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
title('(b)', FontSize = 15)
axis tight;

xlabel('\rho',FontSize=20);
%ylabel('E',FontSize=20);
l = legend('Location', 'northeast',FontSize=15);
title(l,['\beta = ' num2str(beta2_fit)]);
hold off

subplot(1,3,3)
hold on
plot(rho_fit, E_fit(N_b_fit,:),DisplayName = 'Fitted',Color = color3); 
plot(rho_fit, E_num(N_b_fit,:),DisplayName = 'Numerical',LineStyle=":",LineWidth=2,Color = color3); 
plot(rho_ana, E_ana(N_b_ana,:),DisplayName = 'Analytical',LineStyle = '--', LineWidth=2, Color = color3);

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
title('(c)', FontSize = 15)
axis tight;

xlabel('\rho',FontSize=20);
%ylabel('E',FontSize=20);
l = legend('Location', 'northeast',FontSize=15);
title(l,['\beta = ' num2str(beta3_fit)]);
hold off


%Contours for energy
Fig_brEf = figure; 
%sgtitle('Variation of energy vs \rho and \beta', FontSize = 20)

subplot(1,2,1)
contourf(rrho_fit,bb_fit, E_fit, 100, LineStyle = 'None'); 
c = colorbar; 

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
c.FontSize = 12;

xlabel('\rho',FontSize=20);
ylabel('\beta',FontSize=20);
title('(a) Numerical E',FontSize=20);
hold off

subplot(1,2,2)
contourf(rrho_ana,bb_ana, E_ana, 100, LineStyle = 'None'); 
c = colorbar;

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
c.FontSize = 12;

xlabel('\rho', FontSize=20);
ylabel('\beta', FontSize=20);
title('(b) Analytical E',FontSize=20);
hold off

%% Plotting l_x (and l_y)
%Contours
Fig_brl = figure;
subplot(1,2,1);
hold on
contourf(rrho_ana,bb_ana, lx_array, 100, LineStyle = 'None'); 

c = colorbar;
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
c.FontSize = 12;

xlabel('\rho',FontSize=20);
ylabel('\beta',FontSize=20);
title('(a) \lambda_x',FontSize=20);
hold off

subplot(1,2,2)
hold on
contourf(rrho_ana,bb_ana, ly_array, 100, LineStyle = 'None'); 

c = colorbar;
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
c.FontSize = 12;

xlabel('\rho',FontSize=20);
ylabel('\beta',FontSize=20);
title('(b) \lambda_y',FontSize=20);
hold off

%l_x in function of rho
Fig_rlx = figure;
%sgtitle('\lambda_x', FontSize = 20)

subplot(1,2,1)
hold on
plot(rho_ana, lx_array(1,:),DisplayName = ['Num ' num2str(beta1_a)], Color = color1); 
plot(rho_ana, lx_array(ceil(N_b_ana/2),:),DisplayName = ['Num ' num2str(beta2_a)], Color = color2);
plot(rho_ana, lx_array(N_b_ana,:),DisplayName = ['Num ' num2str(beta3_a)], Color = color3); 

plot(rho_ana, lx_array_raw(1,:),'--',DisplayName = ['Ana ' num2str(beta1_a)],LineWidth=2, Color = color1);
plot(rho_ana, lx_array_raw(ceil(N_b_ana/2),:),'--',DisplayName = ['Ana ' num2str(beta2_a)],LineWidth=2, Color = color2);
plot(rho_ana, lx_array_raw(N_b_ana,:),'--',DisplayName = ['Ana ' num2str(beta3_a)],LineWidth=2, Color = color3); 


ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
title('(a)', FontSize = 15)
top = lx_array(1,1)+0.0005;
axis tight;
ylim([-inf top]);

xlabel('\rho',FontSize=20);
ylabel('\lambda_x',FontSize=20);
l = legend('Location', 'southeast',FontSize=15);
title(l,'\beta');
hold off



%l_x in function of beta

subplot(1,2,2)
hold on
plot(beta_ana, lx_array(:,1),DisplayName = ['Num ' num2str(rho1_a)], Color = color1); 
plot(beta_ana, lx_array(:,ceil(N_r_ana/2)),DisplayName = ['Num ' num2str(rho2_a)], Color = color2); 
plot(beta_ana, lx_array(:,N_r_ana),DisplayName = ['Num ' num2str(rho3_a)], Color = color3); 

plot(beta_ana, lx_array_raw(:,1),'--', DisplayName = ['Ana ' num2str(rho1_a)],LineWidth=2, Color = color1); 
plot(beta_ana, lx_array_raw(:,ceil(N_r_ana/2)),'--',DisplayName = ['Ana ' num2str(rho2_a)],LineWidth=2, Color = color2); 
plot(beta_ana, lx_array_raw(:,N_r_ana),'--',DisplayName = ['Ana ' num2str(rho3_a)],LineWidth=2, Color = color3); 

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
title('(b)', FontSize = 15)
axis tight;


xlabel('\beta',FontSize=20);
ylabel('\lambda_x',FontSize=20);
l = legend('Location', 'southwest',FontSize=15);
title(l,'\rho ');
hold off


%l_y in function of rho
Fig_rly = figure;

%sgtitle('\lambda_y', FontSize = 20)

subplot(1,2,1)
hold on
plot(rho_ana, ly_array(1,:),DisplayName = ['Num ' num2str(beta1_a)], Color = color1); 
plot(rho_ana, ly_array(ceil(N_b_ana/2),:),DisplayName = ['Num ' num2str(beta2_a)], Color = color2);
plot(rho_ana, ly_array(N_b_ana,:),DisplayName = ['Num ' num2str(beta3_a)], Color = color3); 

plot(rho_ana, ly_array_raw(1,:),'--',DisplayName = ['Ana ' num2str(beta1_a)],LineWidth=2, Color = color1);
plot(rho_ana, ly_array_raw(ceil(N_b_ana/2),:),'--',DisplayName = ['Ana ' num2str(beta2_a)],LineWidth=2, Color = color2);
plot(rho_ana, ly_array_raw(N_b_ana,:),'--',DisplayName = ['Ana ' num2str(beta3_a)],LineWidth=2, Color = color3); 


ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
title('(a)', FontSize = 15)
top = ly_array(1,1)+0.001;
axis tight;
ylim([-inf top])


xlabel('\rho',FontSize=20);
ylabel('\lambda_y',FontSize=20);
l = legend('Location', 'northwest',FontSize=15);
title(l,'\beta');
hold off

%l_y in function of beta

subplot(1,2,2)
hold on
plot(beta_ana, ly_array(:,1),DisplayName = ['Num ' num2str(rho1_a)], Color = color1); 
plot(beta_ana, ly_array(:,ceil(N_r_ana/2)),DisplayName = ['Num ' num2str(rho2_a)], Color = color2); 
plot(beta_ana, ly_array(:,N_r_ana),DisplayName = ['Num ' num2str(rho3_a)], Color = color3); 

plot(beta_ana, ly_array_raw(:,1),'--', DisplayName = ['Ana ' num2str(rho1_a)],LineWidth=2, Color = color1); 
plot(beta_ana, ly_array_raw(:,ceil(N_r_ana/2)),'--',DisplayName = ['Ana ' num2str(rho2_a)],LineWidth=2, Color = color2); 
plot(beta_ana, ly_array_raw(:,N_r_ana),'--',DisplayName = ['Ana ' num2str(rho3_a)],LineWidth=2, Color = color3); 

ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
title('(b)', FontSize = 15)
axis tight;


xlabel('\beta',FontSize=20);
ylabel('\lambda_y',FontSize=20);
l = legend('Location', 'southwest',FontSize=15);
title(l,'\rho ');
hold off