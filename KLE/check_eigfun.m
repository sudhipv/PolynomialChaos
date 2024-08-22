


%%% Check
clear all
b = 1;
c = 1/b;
sigma_g = 1;
x = -0.5:0.0001:0.5;

eigfun = dlmread('eigfun_check.dat');
g_x = dlmread('g_x_check.dat');
lambda = dlmread('lambda_check.dat');



for k = 1:10

    g_x_cal = zeros(1,length(x));
    
for j = 1:length(x)
    
for i = 1:length(x)

% Simpsons Rule

% % % simp_1   = eigfun(k,i) * exp( -c * abs(x(j)-x(i)));
% % % simp_2  = sigma_g^2 * (4 * eigfun(k,i+1) *exp( -c * abs(x(j)-x(i+1))) ...
% % %     + 2 * eigfun(k,i+2) *exp( -c * abs(x(j)-x(i+2))));
% % % 
% % % if(i ==1 || i == (length(x)-2))
% % %     expo_cov = simp_2 + simp_1 ;
% % % else
% % %     expo_cov = simp_2;
% % % end
% % %     
% % % g_x_cal(:,j) = g_x_cal(:,j) + (0.0001/3)*(expo_cov);

% Rectangular Method
expo_cov = sigma_g^2 * exp( -c * abs(x(j)-x(i)));
g_x_cal(:,j) = g_x_cal(:,j) + (0.0001*eigfun(k,i)*expo_cov);


end

end

figure(k)
plot(x,g_x(k,:))
hold on
plot(x,g_x_cal/sqrt(lambda(k)))
ylim([-inf inf])

end