
% PCE expansion for lognormal random variable using NISP
% using MCS
clear all

mu_g = 1;
sigma_g = .3;
L = 6;

n = 10000;
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


% % % a0 = (sum(Y.*psi_0)/n)/(sum(psi_0.*psi_0)/n)
% % % a1 = (sum(Y.*psi_1)/n)/(sum(psi_1.*psi_1)/n)
% % % a2 = (sum(Y.*psi_2)/n)/(sum(psi_2.*psi_2)/n)
% % % a3 = (sum(Y.*psi_3)/n)/(sum(psi_3.*psi_3)/n)
% % % a4 = (sum(Y.*psi_4)/n)/(sum(psi_4.*psi_4)/n)
% % % a5 = (sum(Y.*psi_5)/n)/(sum(psi_5.*psi_5)/n)

a0 = (sum(Y.*psi_0)/n)/1
a1 = (sum(Y.*psi_1)/n)/1
a2 = (sum(Y.*psi_2)/n)/2
a3 = (sum(Y.*psi_3)/n)/6
a4 = (sum(Y.*psi_4)/n)/24
a5 = (sum(Y.*psi_5)/n)/120


n_MC = 300000;
eta_MC = randn(1,n_MC);

n_MCS = 50000
eta_MCS = randn(1,n_MCS);
Y_MCS = zeros(1,n_MCS);

Y_MCS = 1./exp(mu_g + sigma_g.*eta_MCS);
fprintf('Mean of lognormal random variable by direct MCS')
mean(Y_MCS)

%%% Generating pdf for NISP
psi_0_pdf = ones(1,n_MC);
psi_1_pdf = eta_MC;
psi_2_pdf = (eta_MC.^2 -1);
psi_3_pdf = eta_MC.^3-3.*eta_MC;
psi_4_pdf = eta_MC.^4-6.*eta_MC.^2+3;
psi_5_pdf = eta_MC.^5-10.*eta_MC.^3+15.*eta_MC;
Y_m = a0 .* psi_0_pdf + a1 .* psi_1_pdf +a2 .* psi_2_pdf +a3 .* psi_3_pdf +a4 .* psi_4_pdf +a5 .* psi_5_pdf;

fprintf('mean of PCE log normal')
mean(Y_m)





figure(1)
[fi,xi] = ksdensity(Y_m);
plot(xi,fi)
hold on
[fi_MC,xi_MC] = ksdensity(Y_MCS);
plot(xi_MC,fi_MC)
% % % histogram(Y_m,'Normalization','pdf')
% % % legend('')
% % % hold on
% % % histogram(Y_MC,'Normalization','pdf')
legend('PCE Approximation','Direct MCS')

a = [a0,a1,a2,a3,a4,a5]

%%% Analytical values
    
    
% % % % % %     mu_l = exp(mu_g + 0.5 * (sigma_g^2)); 
% % % % % %     
% % % % % %     psi_analytical = [1, sigma_g,(sigma_g^2)/2,(sigma_g^3)/6,(sigma_g^4)/24,((sigma_g^5)/120)];
% % % % % %     for i = 1:L
% % % % % %     l_analytical(i) = mu_l * psi_analytical(i);
% % % % % %     end
% % % % % % a
% % % % % % l_analytical

