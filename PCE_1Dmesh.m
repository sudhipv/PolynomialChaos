% PCE with MCS for a 1D mesh

x = 1:0.1:10;

L = length(x);


n = 500;

eta = randn(1,n);

Y = zeros(L,n);

for j = 1:n
  
    Y(:,j) = 1 + eta(:,j)*cos((4*pi/10)*x) + (eta(:,j)^2-1)*cos(8*(pi/10)*x);
    
end

psi_0 = ones(1,n);
psi_1 = eta;
psi_2 = (eta.^2 -1);

% ak = <y psi,k>/<psi,k psi,k>

a0 = zeros(L,1);
a1 = zeros(L,1);
a2 = zeros(L,1);


for i = 1:L

a0(i) = (sum(Y(i,:)'.*psi_0(:))/n)/(sum(psi_0(:).*psi_0(:))/n);
a1(i) = (sum(Y(i,:)'.*psi_1(:))/n)/(sum(psi_1(:).*psi_1(:))/n);
a2(i) = (sum(Y(i,:)'.*psi_2(:))/n)/(sum(psi_2(:).*psi_2(:))/n);

end

figure(1)
plot(cos((4*pi/10)*x)); 
hold on; 
plot(a1,'or');
figure(2)
plot(cos(8*(pi/10)*x)); 
hold on; plot(a2,'*b');




