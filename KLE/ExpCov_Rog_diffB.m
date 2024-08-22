
%%% Program to find out the eigen values and eigen functions for
%%% exponential covaraince kernel

%.....Sudhi Sharma P V ...........


% Trying to plot the function to find the approximate roots

% % % clear all
% % % w = 0.01:.0001:0.5;
% % % b = 114;
% % % a = b/2;
% % % val = zeros(1,length(w));
% % % 
% % % 
% % % for i = 1:length(w)
% % %   
% % %     %%%% ODD
% % % %       val(i) = (1/b - w(i)*tan(w(i)*a));
% % %     
% % %     %%% EVEN
% % %      val(i) = (w(i) + (1/b) * tan(w(i)*a));
% % %    
% % %    
% % % %%% Transcedental equation example
% % % %    val(i) = (w(i) * tan(w(i)) - 1);
% % % 
% % %       
% % % end
% % % 
% % % figure(5)
% % % plot(w,zeros(1,length(w)),'-')
% % % hold on
% % % plot(w,val,'b')
% % % ylim([-5,5])
% % % line(w,0)

% Roots from the graphical method for 0-1
%..............................%
% Odd       Even
% .......................
% 1.3065    3.6733
% 6.5846    9.6318
% 12.7234   15.8345
% 18.9550   22.0819
% 25.2120   28.3455
 
%...............................%

% Using fzero to find correct roots

clear all
% Correlation length
b = 1;
a = b/2; % Half the domain length 
mult = 1;

x = 0:0.01:2*a;


n = 5000;
wroot = zeros(1,n);
init = 0.5;
j_o = 1;
check_omega = 0;
sigma_g = 1;
step = 0;

% Odd terms - Equations are swapped in Roger Ghanem Book
for i= 1:n
    
fun =@(w)(1/b - w*tan(w*a));

[wroot(i) fval(i) flag] = fzero(fun,init);

 
        if(flag ~=1)
            init = init +0.1;
            continue
        end

       if(any(abs(check_omega - wroot(i))>0.01,1))
        
           omega_o(j_o) = wroot(i);
           check_omega(j_o) = omega_o(j_o);
           lambda_o(j_o) = sigma_g^2 * 2*b/(1+b^2*omega_o(j_o)^2);
           init = omega_o(j_o);
           
           if(j_o>1)
                  step = omega_o(j_o) - omega_o(j_o-1);
           else
                 step = 2*pi*mult;
%                   step = 0.01;
                  
           end
       
           j_o = j_o+1;
           init = init +step;
           
           if(j_o == 10)
               break
           end
        
       end
       
    
end

lambda_o
omega_o


% Even Terms
init = 0.1;
check_omega = 0;
j_e = 1;
step =0.1;
ne= 50000;


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
             lambda_e(j_e) = sigma_g^2 * 2*b/(1+b^2*omega_e(j_e)^2);
             init = omega_e(j_e);
             %step = 2*pi*mult;
             step = 2*pi*mult;
             j_e = j_e+1;
             
             if(j_e == 10)
               break
             end
             
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
set(gca, 'YScale', 'log','ytick', [0.0001 0.001 0.01 0.1 1]);


%%%% Relative Partial sum of eigen values%%%%
sum = sum(lambda);
relsum = 0;
RelEnergy= zeros(1,length(lambda));

for i = 1:length(lambda)   

relsum = relsum +lambda(i);

RelEnergy(i) = relsum/sum;

end

figure(5)
plot(1:length(lambda),RelEnergy,'-*')
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
  g_x(i,:) = sqrt(lambda(i))*eigfun(i,:);
    
end

% % file = 'eigfun.dat';
% % save(file,'g_x','-ascii')

y=x;
figure(4)
plot(y,eigfun(1:4,:))
legend


%%% Checking eigen functions 

% integral_L = exp(sigma_g * -1 * mod(x1(i)-x2(i))/b) * eigfun()


% % % % % lambda_01 = [1.1648    0.1481    0.0461    0.0216    0.0124    0.0080    0.0056    0.0041    0.0031    0.0025];
% % % % % lambda_1 = [0.7388    0.1380    0.0451    0.0213    0.0123    0.0079    0.0056    0.0041    0.0031    0.0025];
% % % % % lambda_3 = [0.1868    0.0889    0.0382    0.0197    0.0117    0.0077    0.0054    0.0040    0.0031    0.0025];
% % % % % lambda_5 =  [0.0749    0.0520    0.0293    0.0170    0.0107    0.0073    0.0052    0.0039    0.0030    0.0024];
% % % % % lambda_7 = [0.0394    0.0320    0.0217    0.0141    0.0095    0.0067    0.0049    0.0037    0.0029    0.0023];
% % % % % lambda_10 = [0.0197    0.0176    0.0140    0.0104    0.0076    0.0057    0.0044    0.0034    0.0027    0.0022];
% % % % % 
% % % % % figure(5)
% % % % % plot(1:length(lambda_01),lambda_01,'-*')
% % % % % hold on
% % % % % plot(1:length(lambda_1),lambda_1,'-*')
% % % % % plot(1:length(lambda_3),lambda_3,'-*')
% % % % % plot(1:length(lambda_5),lambda_5,'-*')
% % % % % plot(1:length(lambda_7),lambda_7,'-*')
% % % % % plot(1:length(lambda_10),lambda_10,'-*')






