
%%% Program to find out the eigen values and eigen functions for
%%% exponential covaraince kernel

%.....Sudhi Sharma P V ...........

% Using fzero to find correct roots

clear all
delete 'eigfun.dat'
% Correlation length
b = 1;
a = 0.1; % Half the domain length 
mult = 0.5/a;

x = 0:0.01:0.2;


n = 400;
wroot = zeros(1,n);
init = 1;
j_o = 1;
check_omega = 0;
sigma_g = 1;
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
           lambda_o(j_o) = 2*b/(1+b^2*omega_o(j_o)^2);
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
ne= 500;


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
             lambda_e(j_e) = 2*b/(1+b^2*omega_e(j_e)^2);
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


figure(3)
plot(1:length(lambda),lambda,'-*')
ax = gca;
set(gca, 'YScale', 'linear','ytick', [0.0001 0.001 0.01 0.1 1]);

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
  g_x(i,:) = sqrt(lambda(i))*eigfun(i,:);
    
end

file = 'lambda.dat';
save(file,'lambda','-ascii')

file = 'g_x.dat';
save(file,'g_x','-ascii')

file = 'eigfun.dat';
save(file,'eigfun','-ascii')

y=x;
figure(4)
plot(y,eigfun(1:4,:))
legend


%%% Lognormal PCE coefficients

%Gaussian Process with 0 mean
K_n = 4; % Number of terms in KLE
mu = zeros(1,length(x));
g_0 = mu;
g_i = zeros(1,length(x));

for i = 1:K_n-1
    
    g_i = g_i + 0.5 * g_x(i,:).^2;
    
    
end

L_o = exp(g_0 + g_i);

% PCE coefficients of Lognormal Process with inherent gaussian process


L_x(1,:) = L_o;
L_x(2,:) = L_o .* g_x(1,:);
L_x(3,:) = L_o .* g_x(2,:); 
L_x(4,:) = L_o .* g_x(3,:); 
L_x(5,:) = L_o .* g_x(1,:).^2/2; 
L_x(6,:) = L_o .* g_x(1,:).*g_x(2,:);
L_x(7,:) = L_o .* g_x(1,:).*g_x(3,:); 
L_x(8,:) = L_o .* g_x(2,:).^2/2; 
L_x(9,:) = L_o .* g_x(2,:).*g_x(3,:);
L_x(10,:) = L_o .* g_x(3,:).^2/2; 

figure(5)
plot(x,L_x(1:6,:))
legend












