load data_critical_speed.mat
%This file plots the critical speed graph

%Colors
color1 = "#30a8c6";
color2 = "#b53f26";
color3 = "#097969";

%Plotting
scatter(rho_values,v_crit,'filled', DisplayName='Numerical')

hold on
plot(rho_values, sqrt(1-rho_values.^2),color = color2, DisplayName='Analytical')

ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;

xlabel('\rho',FontSize=20);
ylabel('\beta critical',FontSize=20);
leg = legend('Location', 'southwest',FontSize=15);
title(leg,['\beta = ' num2str(beta1_fit)])

hold off