
% PCE Intrusive for lognormal random variable
clear all
% By direct integration
syms si
a = 1;
p = 6; % Number of terms in expansion 
psi = [1,si,(si^2)-1,si^3-3*si,si^4-6*si^2+3,si^5-10*si^3+15*si];
mu_x = 0;
sigma_x = 1;

for i = 1:p
y(i) = (int(exp(mu_x+sigma_x*si)*psi(i)*(1/sqrt(2*pi))*exp(-0.5*(si^2)),si,-inf,inf))/factorial(i-1);
end

double(y) % Coefficients

% By analytical derived formula (Ref: Waad Subber Thesis)

mu_g = mu_x;
sigma_g = sigma_x;

mu_l = exp(mu_g + 0.5 * sigma_g^2);
l_i = zeros(1,p);

for i = 1:p
    l_i(i) = mu_l * (sigma_g^(i-1))/factorial(i-1);
end

l_i % COefficients



% constructing the solution using random numbers
x = 0:.1:5;
n = 10000;

si_m = randn(1,n);
psi_m = [ones(1,n);si_m;(si_m.^2 -1);(si_m.^3 -3.*si_m);si_m.^4-6.*si_m.^2+3;si_m.^5-10.*si_m.^3+15.*si_m];
Y_m = zeros(1,n);

 for i = 1:p
 Y_m = Y_m + double(y(i)).*psi_m(i,:);
 end

fprintf('mean and variance of lognormal distribution')
[M,V] = lognstat(mu_x,sigma_x)
fprintf('mean of Y(PCE)')
mean(Y_m)
fprintf('variance of Y(PCE)')
var(Y_m)
Y_og = lognpdf(x,mu_x,sigma_x);
plot(x,Y_og,'r')
hold on
histogram(Y_m,'Normalization','pdf')
legend('Log Normal PDF-Matlab','Histogram of PCE')

