% PCE - 2nd order with 2 variables for a 1D mesh

x = 1:0.1:10;

L = length(x);


n = 10000;

eta1 = randn(1,n);
eta2 = randn(1,n);

Y = zeros(L,n);

for j = 1:n
  
    Y(:,j) = 1 + eta1(:,j)*cos((4*pi/10)*x) + eta2(:,j)*cos((8*pi/10)*x) + ...
        (eta1(:,j)^2-1)*cos(4*(pi/10)*x) + (eta2(:,j)^2-1)*cos(8*(pi/10)*x)+eta1(:,j)*eta2(:,j)*cos((16*pi/10)*x)+...
        (eta1(:,j)^3-(3*eta1(:,j)))*cos(4*(pi/10)*x)+ (eta2(:,j)^3-(3*eta2(:,j)))*cos(32*(pi/10)*x)+...
    eta2(:,j)*(eta1(:,j)^2-1)*cos(8*(pi/10)*x)+eta1(:,j)*(eta2(:,j)^2-1)*cos(16*(pi/10)*x);
    
end

% for i = 1:L
% sumY(i)= sum(Y(i,:))/n;
% end

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

a0 = zeros(L,1);
a1 = zeros(L,1);
a2 = zeros(L,1);
a3 = zeros(L,1);
a4 = zeros(L,1);
a5 = zeros(L,1);
a6 = zeros(L,1);
a7 = zeros(L,1);
a8 = zeros(L,1);
a9 = zeros(L,1);


for i = 1:L

a0(i) = (sum(Y(i,:)'.*psi_0(:))/n)/(sum(psi_0(:).*psi_0(:))/n);
a1(i) = (sum(Y(i,:)'.*psi_1(:))/n)/(sum(psi_1(:).*psi_1(:))/n);
a2(i) = (sum(Y(i,:)'.*psi_2(:))/n)/(sum(psi_2(:).*psi_2(:))/n);

a3(i) = (sum(Y(i,:)'.*psi_3(:))/n)/(sum(psi_3(:).*psi_3(:))/n);
a4(i) = (sum(Y(i,:)'.*psi_4(:))/n)/(sum(psi_4(:).*psi_4(:))/n);
a5(i) = (sum(Y(i,:)'.*psi_5(:))/n)/(sum(psi_5(:).*psi_5(:))/n);

a6(i) = (sum(Y(i,:)'.*psi_6(:))/n)/(sum(psi_6(:).*psi_6(:))/n);
a7(i) = (sum(Y(i,:)'.*psi_7(:))/n)/(sum(psi_7(:).*psi_7(:))/n);
a8(i) = (sum(Y(i,:)'.*psi_8(:))/n)/(sum(psi_8(:).*psi_8(:))/n);
a9(i) = (sum(Y(i,:)'.*psi_9(:))/n)/(sum(psi_9(:).*psi_9(:))/n);

end

% Recalculating Y
Ycheck = zeros(L,n);
for i = 1:L
Ycheck(i,:) = a0(i) * psi_0 + a1(i) * psi_1 + a2(i) * psi_2 + a3(i) * psi_3 + a4(i) * psi_4 +...
    a5(i) * psi_5 + a6(i) * psi_6 + a7(i) * psi_7 + a8(i) * psi_8 + a9(i) * psi_9;  
end
% 
% for i = 1:L
% sumYcheck(i)= sum(Ycheck(i,:))/n;
% end


figure(1)
plot(cos((4*pi/10)*x)); 
hold on; 
plot(a1,'or');

figure(2)
plot(cos(8*(pi/10)*x)); 
hold on; plot(a2,'*b');

figure(3)
plot(cos(4*(pi/10)*x)); 
hold on; plot(a3,'*b');

figure(4)
plot(cos(16*(pi/10)*x)); 
hold on; plot(a4,'*b');

figure(5)
plot(cos(8*(pi/10)*x)); 
hold on; plot(a5,'*b');

figure(6)
plot(cos(4*(pi/10)*x)); 
hold on; plot(a6,'*b');

figure(7)
plot(cos(8*(pi/10)*x)); 
hold on; plot(a7,'*b');

figure(8)
plot(cos(16*(pi/10)*x)); 
hold on; plot(a8,'*b');

figure(9)
plot(cos(32*(pi/10)*x)); 
hold on; plot(a9,'*b');

% Plotting complete Sloution
% figure(100)
% plot(x,sumY,'r')
% hold on
% plot(x,sumYcheck,'b')

