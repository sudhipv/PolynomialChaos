%% Lognormal Process with inherent Gaussian Process having exponential covaraince kernel using MCS

clear all
% mean and covariance of the 2nd order gaussian process
mu_g = 0;
sigma_g = .1;
L = 6; % Number of terms in expansion

syms si_1
syms si_2

psi_si = [1,si_1,si_2,(si_1^2)-1,si_1*si_2,(si_2^2)-1]
% %     si_1^3-3*si_1,si_2 *((si_1^2)-1),si_1 *((si_2^2)-1),...
% %     si_2^3-3*si_2];

l_0 = exp(mu_g + 0.5 * L * sigma_g^2);

a = zeros(1,L);

a(1) = l_0;
a(2) = l_0 * sigma_g/factorial(1);
a(3) = l_0 * sigma_g/factorial(1);
a(4) = l_0 * sigma_g^2 /factorial(2);
a(5) = l_0 * sigma_g^2 /factorial(1);
a(6) = l_0 * sigma_g^2 /factorial(2);
a(7) = l_0 * sigma_g^3 /factorial(3);
a(8) = l_0 * sigma_g^3 /factorial(2); 
a(9) = l_0 * sigma_g^3 /factorial(2);
a(10) = l_0 * sigma_g^3 /factorial(3);