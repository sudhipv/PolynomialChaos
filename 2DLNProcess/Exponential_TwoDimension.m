

%% Exponential Covariance - Two dimension

% ......Sudhi Sharma P V ..............

%%Taking the tensor Product in Both directions


clear all
delete 'eigfun_2D.dat'
b = 1;
a = 0.5; % Half the domain length 
mult = 0.5/a;

%KLE random dimension
KLE_dim = 3;
% PCE input order
ord_in = 2;
% Number of Spectral Terms in expansion of input random coefficient
n_inpce = factorial(ord_in + KLE_dim)/(factorial(ord_in)*factorial(KLE_dim))


%Underlying Gaussian parameters
mu_g = 0;
sigma_g = .3;

n = 100;
wroot = zeros(1,n);
init = 1;
j_o = 1;
check_omega = 0;

step = 0;

% Odd terms - Equations are swapped in Roger Ghanem Book
for i= 1:n
    
fun =@(w)(1/b - w*tan(w*a));

[wroot(i) fval(i) flag] = fzero(fun,init);

 
        if(flag ~=1)
            init = init +0.01;
            continue
        end

       if(any(abs(check_omega - wroot(i))>0.01,1))
        
           omega_o(j_o) = wroot(i);
           check_omega(j_o) = omega_o(j_o);
           lambda_o(j_o) = sigma_g*2*b/(1+b^2*omega_o(j_o)^2);
           init = omega_o(j_o);
           
           if(j_o>1)
                  step = omega_o(j_o) - omega_o(j_o-1);
           else
                  step = 2*pi*mult;
           end
       
           j_o = j_o+1;
           init = init +step;
        
       end
       
    
end

lambda_o
omega_o


% Even Terms
init = 3.60;
check_omega = 0;
j_e = 1;
step =0;
ne= 100;


for i= 1:ne
    
fun_e =@(w)(w + (1/b)*tan(w*a));


[wroot_e(i) fval_e(i) flag_e] = fzero(fun_e,init);
 
        
% Condition added since froot gives a flag of 1 even after giving singular
% root value
if(wroot_e(i)>1)
    
    if(flag_e ~=1)
            init = init +0.1;
            continue
    end
        
         if(any(abs(check_omega - wroot_e(i))>0.01,1))
             omega_e(j_e) = wroot_e(i);
             check_omega(j_e) = omega_e(j_e);
             lambda_e(j_e) = sigma_g*2*b/(1+b^2*omega_e(j_e)^2);
             init = omega_e(j_e);
             step = 2*pi*mult;
             j_e = j_e+1;
             init = init + step;
         else
            init = init +0.1; 
         end
         
else
    init = init +0.1;
end

    
end

lambda_e  
omega_e

% Ordered Omega and Lambda according to odd and even 
n = 10;
n_e = 1;
n_o = 1;
del_odd = 0;
del_eve = 0;

for i= 1:n
       
   if(mod(i,2) ~=0)
       lambda(i) = lambda_o(n_o);
       omega(i) = omega_o(n_o);
       n_o = n_o + 1;
   
   else
       lambda(i) = lambda_e(n_e);
       omega(i) = omega_e(n_e);
       n_e = n_e + 1;

   end
       
end
    
 lambda  = sort(unique(lambda),'descend')
 omega   = sort(unique(omega))


k = 1;
for i = 1:length(lambda)
    for j= 1:length(lambda)
        
        lambda_n(k,1) = lambda(i)*lambda(j);
        lambda_n(k,2) = i;
        lambda_n(k,3) = j;
        k = k+1;
    end 
end
[B,I] = sort(lambda_n(:,1),'descend');
lambda_multi = lambda_n(I,:,:)


figure(3)
plot(1:length(lambda_multi(1:100,1)),lambda_multi(1:100,1),'*')
ax = gca;
set(gca, 'YScale', 'linear','ytick', [0.0001 0.001 0.01 0.1 1]);

%%%% Relative Partial sum of eigen values%%%%
sum = sum(lambda_multi(1:100,1),1);
relsum = 0;
RelEnergy= zeros(1,100);

for i = 1:100   

relsum = relsum +lambda_multi(i,1);

RelEnergy(i) = relsum/sum;

end

figure(5)
plot(1:100,RelEnergy,'-*')
title('Relative partial sum of eigen values for 2D exponential covariance function')
xlabel('Eigen indices')
ylabel('Relative partial sum')


% Finding the eigen functions using tensor product in the order of
% descending lambda

x = -0.5:0.01:0.5;
y=-0.5:0.01:0.5;
eigfun = zeros(length(omega),length(x));

for i = 1:length(omega)
   
    if(mod(i,2)~=0) % Odd
        den_odd = sqrt(a + sin(2*omega(i)*a)/(2*omega(i)));
        eigfun(i,:) = cos(omega(i).*x)/den_odd;
    else 
        den_even = sqrt(a - sin(2*omega(i)*a)/(2*omega(i)));
        eigfun(i,:) = sin(omega(i).*x)/den_even;
    end
    
    eigfun(i,:) = eigfun(i,:);

    
end

% 2D Eigen Functions and KLE coefficient terms

file = 'eigfun_2D.dat';
for i = 1:n_inpce
eigfun_multi(:,:,i) = eigfun(lambda_multi(i,2),:)'*eigfun(lambda_multi(i,3),:);
g_x_2D(:,:,i) = sqrt(lambda_multi(i,1))*eigfun(lambda_multi(i,2),:)'*eigfun(lambda_multi(i,3),:);
dlmwrite(file,g_x_2D,'-append')
g_x_2D(:,:,i) = g_x_2D(:,:,i).*sigma_g;

end



for i=1:n_inpce
figure(i)
surf(x,y,eigfun_multi(:,:,i))
xlabel('X')
ylabel('Y')
zlabel('\surd(\lambda_i_x * \lambda_j_y) g_i(x) g_j(y)')
end



%%% Lognormal PCE coefficients 2Dimension

%Gaussian Process with 0 mean
% % mu = ones(length(x),length(y));
% % g_0 = mu_g*mu;
% % g_i = zeros(length(x),length(y));
% % 
% % for i = 1:KLE_dim
% %     
% %     g_i = g_i + 0.5 * g_x_2D(:,:,i).^2;
% %     
% %     
% % end

%Finding Mean from Number of terms used
% % % L_o = exp(g_0 + g_i)


%   Finding Mean of Log Normal Expansion directly

L_o = exp(mu_g + 0.5*sigma_g^2);

% PCE coefficients of Lognormal Process with inherent gaussian process


L_x_2D(:,:,1) = L_o.*ones(length(x),length(y)); % 1
L_x_2D(:,:,2) = L_o .* g_x_2D(:,:,1); % xi1
L_x_2D(:,:,3) = L_o .* g_x_2D(:,:,2); % xi2
L_x_2D(:,:,4) = L_o .* g_x_2D(:,:,3); % xi3
L_x_2D(:,:,5) = L_o .* g_x_2D(:,:,1).^2/2; % xi1^2-1
L_x_2D(:,:,6) = L_o .* g_x_2D(:,:,1).*g_x_2D(:,:,2); % xi1*xi2;
L_x_2D(:,:,7) = L_o .* g_x_2D(:,:,1).*g_x_2D(:,:,3); % xi1*xi3; 
L_x_2D(:,:,8) = L_o .* g_x_2D(:,:,2).^2/2; % xi2^2-1;
L_x_2D(:,:,9) = L_o .* g_x_2D(:,:,2).*g_x_2D(:,:,3); % xi2*xi3; 
L_x_2D(:,:,10) = L_o .* g_x_2D(:,:,3).^2/2; % xi3^2-1;
 

for i = 1:n_inpce
figure(10+i)
surf(x,y,L_x_2D(:,:,i))
end
legend



