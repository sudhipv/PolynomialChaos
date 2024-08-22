
%%%%%%%%% Non Intrusive %%%%%%%%%

% random dimension of 2

mu_g = 0;
sigma_g(1) = 1;
sigma_g(2) = 1;

n = 100000;
eta1 = randn(1,n);
eta2 = randn(1,n);
L = 10;

Y = zeros(1,n);
Y = exp(mu_g + (sigma_g(1).*eta1) + (sigma_g(2).*eta2));
fprintf('Mean of log normal RV by direct MCS');
mean(Y)

psi_0 = ones(1,n);
psi_1 = eta1;
psi_2 = eta2;
psi_3 = (eta1.^2 -1);
psi_4 = (eta1.*eta2);
psi_5 = (eta2.^2 -1);
psi_6 = (eta1.^3 -(3.*eta1));
psi_7 = eta2.*(eta1.^2 -1);
psi_8 = eta1.*(eta2.^2 -1);
psi_9 = (eta2.^3 -(3.*eta2));
% ak = <y psi,k>/<psi,k psi,k>

a = zeros(1,L);

a(1) = (sum(Y.*psi_0)/n)/(sum(psi_0.*psi_0)/n);
a(2) = (sum(Y.*psi_1)/n)/(sum(psi_1.*psi_1)/n);
a(3) = (sum(Y.*psi_2)/n)/(sum(psi_2.*psi_2)/n);

a(4) = (sum(Y.*psi_3)/n)/(sum(psi_3.*psi_3)/n);
a(5) = (sum(Y.*psi_4)/n)/(sum(psi_4.*psi_4)/n);
a(6) = (sum(Y.*psi_5)/n)/(sum(psi_5.*psi_5)/n);

a(7) = (sum(Y.*psi_6)/n)/(sum(psi_6.*psi_6)/n);
a(8) = (sum(Y.*psi_7)/n)/(sum(psi_7.*psi_7)/n);
a(9) = (sum(Y.*psi_8)/n)/(sum(psi_8.*psi_8)/n);
a(10) = (sum(Y.*psi_9)/n)/(sum(psi_9.*psi_9)/n);


%%% Analytical values
    
    
    mu_l = exp(mu_g + 0.5 * (sigma_g(1)^2+sigma_g(2)^2)); 
    
    psi_analytical = [1, sigma_g(1),sigma_g(2),(sigma_g(1)^2)/2, sigma_g(1)*sigma_g(2), (sigma_g(2)^2)/2, ...
       (sigma_g(1)^3)/6, (sigma_g(1)^2*sigma_g(2))/2,(sigma_g(1)*sigma_g(2)^2)/2,(sigma_g(2)^3)/6];
    for i = 1:L
    l_analytical(i) = mu_l * psi_analytical(i);
    end
a
l_analytical

