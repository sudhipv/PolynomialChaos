

%%%% Non intrusive Method %%%%%%
% Log Normal Process with Exponential Covariance 1D 

clear all


x = -0.5:0.01:0.5;

% eigfun = sqrt(lambda)*phi(x)
% This is generated from example "ExpoentialCovariance_Roger",
% So any changes pertaining to eigenfunctions and eigen values must be
% altered in that example first to get matching results.
eigfun = load('g_x.dat'); 


%Gaussian Process with 0 mean
P = 4; % Number of terms in KLE
n = 10000;
L = zeros(length(x),n);

mu = zeros(1,length(x));
g_0 = mu;



 
% Generating P random variables of size (1,n)
for i = 1:P-1
xi(i,:) = randn(1,n);
end

% Generating 'n' realisations of the Log normal random process
for j = 1:n
   
    g_x = zeros(1,length(x));
    
for i = 1:P-1
    
    g_x = g_x + xi(i,j).*eigfun(i,:);
    
end

L(:,j) = exp(g_0 + g_x);

end


% Ensemble Mean of Lognormal Process through Direct MCS
for i = 1:length(x)
   L_mean(i) = sum(L(i,:))/n;
end

figure(1)
plot(x,L_mean,'k*-.')
hold on


% PCE Coefficients of Lognormal Process through MCS

psi_0 = ones(1,n);
psi_1 = xi(1,:);
psi_2 = xi(2,:);
psi_3 = xi(3,:);
psi_4 = xi(1,:).^2-1;
psi_5 = xi(1,:).*xi(2,:);
psi_6 = xi(1,:).*xi(3,:);
psi_7 = xi(2,:).^2-1;
psi_8 = xi(2,:).*xi(3,:);
psi_9 = xi(3,:).^2-1;

% % % psi_4 = xi(4,:);
% % % psi_5 = xi(5,:);
% % % psi_6 = xi(1,:).^2-1;
% % % psi_7 = xi(2,:).^2-1;
% % % psi_8 = xi(3,:).^2-1;
% % % psi_9 = xi(4,:).^2-1;
% % % psi_10 = xi(5,:).^2-1;
% % % psi_11 = xi(1,:).*xi(2,:);
% % % psi_12 = xi(1,:).*xi(3,:);
% % % psi_13 = xi(1,:).*xi(4,:);
% % % psi_14 = xi(1,:).*xi(5,:);
% % % psi_15 = xi(2,:).*xi(3,:);


for i = 1:length(x)

a(1,i) = (sum(L(i,:).*psi_0)/n)/(sum(psi_0.*psi_0)/n);
a(2,i) = (sum(L(i,:).*psi_1)/n)/(sum(psi_1.*psi_1)/n);
a(3,i) = (sum(L(i,:).*psi_2)/n)/(sum(psi_2.*psi_2)/n);
a(4,i) = (sum(L(i,:).*psi_3)/n)/(sum(psi_3.*psi_3)/n);
a(5,i) = (sum(L(i,:).*psi_4)/n)/(sum(psi_4.*psi_4)/n);
a(6,i) = (sum(L(i,:).*psi_5)/n)/(sum(psi_5.*psi_5)/n);
a(7,i) = (sum(L(i,:).*psi_6)/n)/(sum(psi_6.*psi_6)/n);
a(8,i) = (sum(L(i,:).*psi_7)/n)/(sum(psi_7.*psi_7)/n);
a(9,i) = (sum(L(i,:).*psi_8)/n)/(sum(psi_8.*psi_8)/n);
a(10,i) = (sum(L(i,:).*psi_9)/n)/(sum(psi_9.*psi_9)/n);


% % % a(11,i) = (sum(L(i,:).*psi_10)/n)/(sum(psi_10.*psi_10)/n);
% % % a(12,i) = (sum(L(i,:).*psi_11)/n)/(sum(psi_11.*psi_11)/n);
% % % a(13,i) = (sum(L(i,:).*psi_12)/n)/(sum(psi_12.*psi_12)/n);
% % % a(14,i) = (sum(L(i,:).*psi_13)/n)/(sum(psi_13.*psi_13)/n);
% % % a(15,i) = (sum(L(i,:).*psi_14)/n)/(sum(psi_14.*psi_14)/n);
% % % a(16,i) = (sum(L(i,:).*psi_15)/n)/(sum(psi_15.*psi_15)/n);


end

for i = 1:6
plot(x,a(i,:))
hold on
end
legend

%%%% Marginal PDF through MCS

L_pce = a(1,:)'.*psi_0 +a(2,:)'.*psi_1 + a(3,:)'.*psi_2 + a(4,:)'.*psi_3 +a(5,:)'.*psi_4 +a(6,:)'.*psi_5...
    + a(7,:)'.*psi_6 + a(8,:)'.*psi_7 + a(9,:)'.*psi_8 + a(10,:)'.*psi_9; 

% % % + a(11,:)'.*psi_10 + a(12,:)'.*psi_11 +...
% % % a(13,:)'.*psi_12 +a(14,:)'.*psi_13 + a(15,:)'.*psi_14 +  a(16,:)'.*psi_15 ; 

% % % % At x = 0;
% % % x1 = ceil(length(x)/2);
% % % figure(2)
% % % histogram(L_pce(x1,:),'Normalization','pdf')
% % % hold on
% % % histogram(L(x1,:),'Normalization','pdf')
% % % legend('Marginal PDF from PCE','Marginal PDF from MCS')
% % % 
% % % 
% % % % At x = -0.25;
% % % figure(3)
% % % x2 = ceil(length(x)/4);
% % % histogram(L_pce(x2,:),'Normalization','pdf')
% % % hold on
% % % histogram(L(x2,:),'Normalization','pdf')
% % % legend('Marginal PDF from PCE','Marginal PDF from MCS')
% % % 
% % % % At x = 0.25;
% % % figure(4)
% % % x3 = ceil(length(x)*3/4);
% % % histogram(L_pce(x3,:),'Normalization','pdf')
% % % hold on
% % % histogram(L(x3,:),'Normalization','pdf')
% % % legend('Marginal PDF from PCE','Marginal PDF from MCS')

% Kernel Density estimates
x1 = ceil(length(x)/2);
x2 = ceil(length(x)/4);
x3 = ceil(length(x)*3/4);
[f_0,xi_0] = ksdensity(L_pce(x1,:));
[f_mcs,xi_mcs] = ksdensity(L(x1,:));
figure(6)
plot(xi_0,f_0,xi_mcs,f_mcs)
legend('Marginal PDF from PCE','Marginal PDF from MCS')

[f_L,xi_L] = ksdensity(L_pce(x2,:));
[f_mcs_L,xi_mcs_L] = ksdensity(L(x2,:));
figure(7)
plot(xi_L,f_L,xi_mcs_L,f_mcs_L)
legend('Marginal PDF from PCE','Marginal PDF from MCS')

[f_R,xi_R] = ksdensity(L_pce(x3,:));
[f_mcs_R,xi_mcs_R] = ksdensity(L(x3,:));
figure(8)
plot(xi_R,f_R,xi_mcs_R,f_mcs_R)
legend('Marginal PDF from PCE','Marginal PDF from MCS')
