
clear all
b = 114;
a = b/2;
x = 0:0.1:2*a;

%%%% calculate the omegas by graphical method

%%%%%%% 144 domain size and same correlation length
omega_odd = [0.166, 0.221, 0.2761, 0.3312, 0.386,0.4413];
omega_even = [0.032,0.0845,0.139,0.194];


%%%%%%% 200 domain size and same correlation length
% omega_odd = [0.033,0.063,0.0948,0.126,0.1573];
% omega_even = [0.0184,0.0481,0.0792,0.11042];

for i = 1:size(omega_odd,2)
l_odd(i,1) = 2*b/(1+b^2*omega_odd(i)^2);
end

for i = 1:size(omega_even,2)
l_even(i,1) = 2*b/(1+b^2*omega_even(i)^2);
end



n = 8;
n_e = 1;
n_o = 1;
del_odd = 0;
del_eve = 0;


for i= 1:8
       
   if(mod(i,2) ~=0)
       lam(i) = l_odd(n_o);
       omega(i) = omega_odd(n_o);
       n_o = n_o + 1;
   
   else
       lam(i) = l_even(n_e);
       omega(i) = omega_even(n_e);
       n_e = n_e + 1;

   end
       
end
    
 lam  = sort(lam,'descend')
 omega   = sort(omega)
 
 
figure(3)
plot(1:length(lam),lam,'-*')
ax = gca;
set(gca, 'YScale', 'log','ytick', [0.0001 0.001 0.01 0.1 1]);


%%%% Relative Partial sum of eigen values%%%%
sum = sum(lam);
relsum = 0;
RelEnergy= zeros(1,length(lam));

for i = 1:length(lam)   

relsum = relsum +lam(i);

RelEnergy(i) = relsum/sum;

end

figure(5)
plot(1:length(lam),RelEnergy,'-*')
title('Relative partial sum of eigen values for 1D exponential covariance function')
xlabel('Eigen indices')
ylabel('Relative partial sum')



% Finding the eigen functions
eigfun = zeros(length(omega),length(x));

for i = 1:length(omega)
   
    if(mod(i,2)~=0)
        den_odd = sqrt(a + sin(2*omega(i)*a)/(2*omega(i)));
        eigfun(i,:) = cos(omega(i).*(x-a))/den_odd;
    else 
        den_even = sqrt(a - sin(2*omega(i)*a)/(2*omega(i)));
        eigfun(i,:) = sin(omega(i).*(x-a))/den_even;
    end
    
    eigfun(i,:) = eigfun(i,:);
%     .* sqrt(lambda(i));
  g_x(i,:) = sqrt(lam(i))*eigfun(i,:);
    
end

% % file = 'eigfun.dat';
% % save(file,'g_x','-ascii')

y=x;
figure(4)
plot(y,eigfun(1:4,:))
legend
 
 


%%%%%%%%%%% Calculating Product of 3 dimensional Lambdas . %%%%%%
%%

lambda_200 = load('lambda_200.mat');
lambda_114 = load('lambda_114.mat');

lambda_200 = lambda_200.lambda_200;
lambda_114 = lambda_114.lambda_114;

n = 1;
for i = 1:length(lambda_200)
    for j= 1:length(lambda_114)
        
        for k= 1:length(lambda_200)
        
            lambda_n(n,1) = lambda_200(i)*lambda_114(j)*lambda_200(k);
            lambda_n(n,2) = i;
            lambda_n(n,3) = j;
            lambda_n(n,4) = k;
            n = n+1;
        
        end
    end 
end
[B,I] = sort(lambda_n(:,1),'descend');
lambda_multi = lambda_n(I,:,:)


figure(3)
plot(1:length(lambda_multi(1:100,1)),lambda_multi(1:100,1),'*')
ax = gca;
set(gca, 'YScale', 'linear','ytick', [0.0001 0.001 0.01 0.1 1]);

%%%% Relative Partial sum of eigen values%%%%
sumlamda = sum(lambda_multi(1:100,1),1);
relsum = 0;
RelEnergy= zeros(1,100);

for i = 1:100   

relsum = relsum +lambda_multi(i,1);

RelEnergy(i) = relsum/sumlamda;

end

figure(5)
plot(1:100,RelEnergy,'-*')
title('Relative partial sum of eigen values for 2D exponential covariance function')
xlabel('Eigen indices')
ylabel('Relative partial sum')



%%

a_114 = 57;

a_200 = 100;

lambda_200 = load('lambda_200.mat');
lambda_114 = load('lambda_114.mat');

lambda_200 = lambda_200.lambda_200;
lambda_114 = lambda_114.lambda_114;


omega_200 = load('omega_200.mat');
omega_114 = load('omega_114.mat');

omega_200 = omega_200.omega_200;
omega_114 = omega_114.omega_114;


mult_114 = zeros(length(omega_114),1);

for i = 1:length(omega_114)
   
    if(mod(i,2)~=0)
        mult_114(i) = sqrt(lambda_114(i))/ ( sqrt(a_114 + ( sin(2*omega_114(i)*a_114)/(2*omega_114(i)) ) ) );

    else 
        mult_114(i) = sqrt(lambda_114(i))/( sqrt(a_114 - ( sin(2*omega_114(i)*a_114)/(2*omega_114(i)) ) ) );
    end
    
  
 
end



mult_200 = zeros(length(omega_200),1);

for i = 1:length(omega_200)
   
    if(mod(i,2)~=0)
        mult_200(i) = sqrt(lambda_200(i))/( sqrt(a_200 + ( sin(2*omega_200(i)*a_200)/(2*omega_200(i)) ) ) );

    else 
        mult_200(i) = sqrt(lambda_200(i))/( sqrt(a_200 -  ( sin(2*omega_200(i)*a_200)/(2*omega_200(i)) ) ) );
    end
    
  
 
end












 
 