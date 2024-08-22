
% PCE expansion for a log Normal random variable expanded using two
% gaussian variables

%%%%%%%%%%Intrusive with MCS %%%%%%%%%%%%

clear all
% mean and variance of the 2nd order gaussian process
mu_g = 0;
sigma_g = 1;
L = 10; % Number of terms in expansion

syms si_1
syms si_2

psi_si = [1,si_1,si_2,(si_1^2)-1,si_1*si_2,(si_2^2)-1,si_1^3-3*si_1,si_2 *((si_1^2)-1),si_1 *((si_2^2)-1),...
    si_2^3-3*si_2];

l_0 = exp(mu_g + 0.5 * sigma_g^2);

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

a

l_xtheta =l_0 * (1 * psi_si(1)+ sigma_g * psi_si(2)/factorial(1)+ sigma_g * psi_si(3)/factorial(1)...
    + sigma_g^2 * psi_si(4)/factorial(2) + sigma_g^2 * psi_si(5)/factorial(1) + sigma_g^2 * psi_si(6)/factorial(2)...
    + sigma_g^3 * psi_si(7)/factorial(3) + sigma_g^3 * psi_si(8)/factorial(2) + sigma_g^3 * psi_si(9)/factorial(2)...
    + sigma_g^3 * psi_si(10)/factorial(3)); % For multi dimensions change sigma_g


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 1000;
eta1 = randn(1,n);
eta2 = randn(1,n);
l_theta = double(subs(l_xtheta,{si_1,si_2},{eta1,eta2}));

mean(l_theta)
var(l_theta)
x = 0:.1:10;
Y_og = lognpdf(x,mu_g,sigma_g);
[m,v] = lognstat(mu_g,sigma_g)
plot(x,Y_og,'r')
hold on
histogram(l_theta,'Normalization','pdf')
legend('Log Normal PDF-Matlab','Log Normal process Marginal pdf using PCE and MCS')
