
% PCE expansion for lognormal random variable using NISP
% using MCS
clearvars

for n_l = 1:20

mu_g = 1;
sigma_g = .3;

n = 50000;
eta = randn(1,n);

Y = zeros(1,n);

Y = 1./exp(mu_g + sigma_g.*eta);
% fprintf('Mean of lognormal random variable by direct MCS')
% mean(Y)

% PCE expansion of Y

psi_0 = ones(1,n);
psi_1 = eta;
psi_2 = (eta.^2 -1);
psi_3 = eta.^3-3.*eta;
psi_4 = eta.^4-6.*eta.^2+3;
psi_5 = eta.^5-10.*eta.^3+15.*eta;
psi_6 = eta.^6-15*eta.^4+45*eta.^2-15;
psi_7 = eta.^7-21*eta.^5+105*eta.^3-105*eta;
psi_8 = eta.^8-28*eta.^6+210*eta.^4-420*eta.^2+105;
psi_9 = eta.^9-36*eta.^7+378*eta.^5-1260*eta.^3+945*eta;
psi_10 = eta.^10-45*eta.^8+630*eta.^6-3150*eta.^4+4725*eta.^2-945;

% % % a0 = (sum(Y.*psi_0)/n)/(sum(psi_0.*psi_0)/n)
% % % a1 = (sum(Y.*psi_1)/n)/(sum(psi_1.*psi_1)/n)
% % % a2 = (sum(Y.*psi_2)/n)/(sum(psi_2.*psi_2)/n)
% % % a3 = (sum(Y.*psi_3)/n)/(sum(psi_3.*psi_3)/n)
% % % a4 = (sum(Y.*psi_4)/n)/(sum(psi_4.*psi_4)/n)
% % % a5 = (sum(Y.*psi_5)/n)/(sum(psi_5.*psi_5)/n)

a0 = (sum(Y.*psi_0)/n)/1;
a1 = (sum(Y.*psi_1)/n)/1;
a2 = (sum(Y.*psi_2)/n)/2;
a3 = (sum(Y.*psi_3)/n)/6;
a4 = (sum(Y.*psi_4)/n)/24;
a5 = (sum(Y.*psi_5)/n)/120;
a6 = (sum(Y.*psi_6)/n)/720;
a7 = (sum(Y.*psi_7)/n)/5040;
a8 = (sum(Y.*psi_8)/n)/40320;
a9 = (sum(Y.*psi_9)/n)/362880;
a10 = (sum(Y.*psi_10)/n)/3628800;

n_MC = 300000;
eta_MC = randn(1,n_MC);

%%% Generating pdf for NISP
psi_0_pdf = ones(1,n_MC);
psi_1_pdf = eta_MC;
psi_2_pdf = (eta_MC.^2 -1);
psi_3_pdf = eta_MC.^3-3.*eta_MC;
psi_4_pdf = eta_MC.^4-6.*eta_MC.^2+3;
psi_5_pdf = eta_MC.^5-10.*eta_MC.^3+15.*eta_MC;
psi_6_pdf = eta_MC.^6-15*eta_MC.^4+45*eta_MC.^2-15;
psi_7_pdf = eta_MC.^7-21*eta_MC.^5+105*eta_MC.^3-105*eta_MC;
psi_8_pdf = eta_MC.^8-28*eta_MC.^6+210*eta_MC.^4-420*eta_MC.^2+105;
psi_9_pdf = eta_MC.^9-36*eta_MC.^7+378*eta_MC.^5-1260*eta_MC.^3+945*eta_MC;
psi_10_pdf = eta_MC.^10-45*eta_MC.^8+630*eta_MC.^6-3150*eta_MC.^4+4725*eta_MC.^2-945;

Y_m = a0 .* psi_0_pdf + a1 .* psi_1_pdf +a2 .* psi_2_pdf +a3 .* psi_3_pdf +a4 .* psi_4_pdf +a5 .* psi_5_pdf+...
         a6 .* psi_6_pdf +a7 .* psi_7_pdf+a8 .* psi_8_pdf+a9 .* psi_9_pdf+a10 .* psi_10_pdf;
     
fprintf('mean of PCE log normal')
mean(Y_m)

load MCS_pcvariation.mat
figure(1)
[fi,xi] = ksdensity(Y_m);
plot(xi_MC,fi_MC,xi,fi)

% % % histogram(Y_m,'Normalization','pdf')
% % % legend('')
% % % hold on
% % % histogram(Y_MC,'Normalization','pdf')
legend('Direct MCS','PCE Approximation')

a = [a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10]

PCvar_coeff50000(n_l,:) = a;


end

%%% Analytical values
    
    
% % % % % %     mu_l = exp(mu_g + 0.5 * (sigma_g^2)); 
% % % % % %     
% % % % % %     psi_analytical = [1, sigma_g,(sigma_g^2)/2,(sigma_g^3)/6,(sigma_g^4)/24,((sigma_g^5)/120)];
% % % % % %     for i = 1:L
% % % % % %     l_analytical(i) = mu_l * psi_analytical(i);
% % % % % %     end
% % % % % % a
% % % % % % l_analytical

