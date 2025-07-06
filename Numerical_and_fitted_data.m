%This file should be run before Plotter.m
%This file allows to obtain the data relative
% to the relaxed numerical skyrmion and the fitted model.

%Loading and normalizing the numerical skyrmion
load skyrmion5_61.mat;
[Sx,Sy,Sz] = normalizeSpins(Sx,Sy,Sz);
[Ny,Nx] = size(Sx); %size, must MAKE Nx and Ny ODD
close all;

%Parameters (for beta and rho ranges)
N_beta = 11;
N_rho = 11;
beta_min = 0.01;
beta_max = 0.26;
rho_min = 0.9;
rho_max = 0.98;

%Ranges for rho and beta
beta_values = linspace(beta_min,beta_max,N_beta);
rho_values = linspace(rho_min,rho_max,N_rho);

[rrho,bb] = meshgrid(rho_values,beta_values);

%Basic interaction coefficients
J = 1;
K = 0.05;

%Relaxation time and damping constant
time = 800;
time_step = 200;
aRel = 1;

%Allocating space
E_array = zeros(size(rrho));
E_mod_array = zeros(size(rrho));
R_array = zeros(size(rrho));
w_array = zeros(size(rrho));
deltaE_array = zeros(size(rrho));
%lx_array_num = zeros(size(rrho)); %In case computing lx and ly by
%fitting
%ly_array_num = zeros(size(rrho));

%Computing radius and width for fitted skyrmions, and energy for both
%fitted and numeric skyrmions

for i=1:N_beta
    for j=1:N_rho
        [E_array(i,j),E_mod_array(i,j),deltaE_array(i,j),R_array(i,j),w_array(i,j)] = ...
            DeltaE_Rw(Sx,Sy,Sz,J,K,bb(i,j),rrho(i,j),time,time_step,aRel);
    end
end

%Renaming of values for further purpose (Plotter.m)
E_num = E_array;
E_fit = E_mod_array;
R_fit = R_array;
w_fit = w_array;
deltaE_fit = deltaE_array;

rho_fit = rho_values;
beta_fit = beta_values;
rrho_fit = rrho;
bb_fit = bb;
N_r_fit = N_rho;
N_b_fit = N_beta;



%Useful values of beta and rho for further use
rho1_fit = rho_fit(1);
rho2_fit = rho_fit(ceil(N_r_fit/2));
rho3_fit = rho_fit(N_r_fit);

beta1_fit = beta_fit(1);
beta2_fit = beta_fit(ceil(N_b_fit/2));
beta3_fit = beta_fit(N_b_fit);

rho1_fit = round(rho1_fit,2);
rho2_fit = round(rho2_fit,2);
rho3_fit = round(rho3_fit,2);

beta1_fit = round(beta1_fit,2);
beta2_fit = round(beta2_fit,2);
beta3_fit = round(beta3_fit,2);

