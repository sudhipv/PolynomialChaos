clear all
% eigen values and eigen functions for exponential co-variance kernel
x = 0:.001:1;

sigma = 1;
b = 1;

for i = 1:length(x)
    for j = 1:length(x)
K(i,j) = sigma^2 * exp(-1*abs((x(i)-x(j)))/b);
    end
end

[V,D] = eig(K);
e = diag(D);
lambda_sort = sort(diag(D),'descend');
lambda_norm =lambda_sort/sum(lambda_sort);
figure(1)
plot(1:10,lambda_norm(1:10),'*')
set(gca, 'YScale', 'log')

figure(2)
for i = length(x):-1:1
plot(x,sqrt(e(i)).*V(:,i))
hold on
end
legend



% %%%% Relative Partial sum of eigen values%%%%
% sum = sum(lambda);
% relsum = 0;
% RelEnergy= zeros(1,length(lambda));
% 
% for i = 1:length(lambda)   
% 
% relsum = relsum +lambda(i);
% 
% RelEnergy(i) = relsum/sum;
% 
% end
% 
% figure(5)
% plot(1:length(lambda),RelEnergy,'-*')
% title('Relative partial sum of eigen values for 1D exponential covariance function')
% xlabel('Eigen indices')
% ylabel('Relative partial sum')
% 






