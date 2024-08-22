
% PCE coefficients through MCS

% Y = X^2

n = 10000;

eta = randn(1,n);

mu = 0;

sigma = 1;

x = mu + sigma * eta;

Y = x.^2;

psi_0 = ones(1,n);
psi_1 = eta;
psi_2 = (eta.^2 -1);

% ak = <y psi,k>/<psi,k psi,k>

a0 = (sum(Y(:).*psi_0(:))/n)/(sum(psi_0(:).*psi_0(:))/n)
a1 = (sum(Y(:).*psi_1(:))/n)/(sum(psi_1(:).*psi_1(:))/n)
a2 = (sum(Y(:).*psi_2(:))/n)/(sum(psi_2(:).*psi_2(:))/n)








