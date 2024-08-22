
% Monte Carlo Simulation and KDE- Kernel Density Estimation

% Y = aX^2 , X ~ N(10,4)

a = 2;
mu = 10;
sigma = 2;

n = 10000;

X = mu + sigma * randn(1,n);
Y = a*X.^2;
y = 0:700;

analytical = (1./(2*a*sqrt(y/a)).* (1/(sqrt(2*pi)*sigma)).*exp(-1*(sqrt(y/a)-mu).^2/(2*sigma^2))) +... 
(1./(2*a*sqrt(y/a)).* (1/(sqrt(2*pi)*sigma)).*exp(-1*(-sqrt(y/a)-mu).^2/(2*sigma^2)));

[f,xi] = ksdensity(Y);

plot(xi,f);
hold on
plot(y,analytical,'k');



