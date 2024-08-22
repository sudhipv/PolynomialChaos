
% PCE - 2nd order, 2 variables for a 2D mesh

clear all

x = 1:0.1:10;
y = 1:0.1:10;

L = length(x);


n = 500;

eta1 = randn(1,n);
eta2 = randn(1,n);

Y = zeros(L,L,n);

for j = 1:n
  
    Y(:,:,j) = 1 + eta1(:,j)*cos((4*pi/10)*x)'*cos((4*pi/10)*y) + eta2(:,j)*cos((8*pi/10)*x)'*cos((8*pi/10)*y) + ...
        (eta1(:,j)^2-1)*cos(4*(pi/10)*x)'*cos(4*(pi/10)*y) + (eta2(:,j)^2-1)*cos(8*(pi/10)*x)'*cos(8*(pi/10)*y)+...
        eta1(:,j)*eta2(:,j)*cos((16*pi/10)*x)'*cos((16*pi/10)*y);
    
end

% 2nd Order

psi_0 = ones(1,n);
psi_1 = eta1;
psi_2 = eta2;
psi_3 = (eta1.^2 -1);
psi_4 = (eta1.*eta2);
psi_5 = (eta2.^2 -1);

a0 = zeros(L,L);
a1 = zeros(L,L);
a2 = zeros(L,L);
a3 = zeros(L,L);
a4 = zeros(L,L);
a5 = zeros(L,L);

for i = 1:L
    
    for j = 1:L

            temp(:,1) = Y(i,j,:);
        
            a0(i,j) = (sum(temp.*psi_0(:))/n)/(sum(psi_0(:).*psi_0(:))/n);
            a1(i,j) = (sum(temp.*psi_1(:))/n)/(sum(psi_1(:).*psi_1(:))/n);
            a2(i,j) = (sum(temp.*psi_2(:))/n)/(sum(psi_2(:).*psi_2(:))/n);

            a3(i,j) = (sum(temp.*psi_3(:))/n)/(sum(psi_3(:).*psi_3(:))/n);
            a4(i,j) = (sum(temp.*psi_4(:))/n)/(sum(psi_4(:).*psi_4(:))/n);
            a5(i,j) = (sum(temp.*psi_5(:))/n)/(sum(psi_5(:).*psi_5(:))/n);
    end
    
end

figure(1)
surf(x,y,a1)
figure(2)
surf(x,y,cos((4*pi/10)*x)'*cos((4*pi/10)*y))
