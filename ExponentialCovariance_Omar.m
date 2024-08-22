% Using Anlytical solutions for exponential co-variance

% TRying to plot the function to find the approximate roots

% % % clear all
% % % format shortG
% % % w = 0.0:.01:20;
% % % b = 1;
% % % val = zeros(1,length(w));
% % % 
% % % 
% % % for i = 1:length(w)
% % %   
% % %     val(i) = (1-b*w(i)*tan(w(i)/2))*(b*w(i)+tan(w(i)/2));
% % %       
% % % end
% % % 
% % % figure(5)
% % % plot(w,val)

% Using fzero to find correct roots

%%% Equations are not generalized for any domain and correlation lengths
%%% only for (0,1) with correlation length 1

clear all
b = 1;
a = b/2;
n = 10;
wroot = zeros(1,n);
init = 1;
j = 1;
check_omega = 0;
step = 0.1;

for i= 1:n
    
% fun =@(w)((1/b - w*tan(a*w))*(w + (1/b)*tan(a*w)));
fun =@(w)((1 - b*w*tan(w/2))*(b*w + tan(w/2)));

init = init + step;
[wroot(i) fval(i) flag] = fzero(fun,init);
 


if(flag ~=1)
    continue
end
    if(abs(wroot(i) - check_omega) > 0.001 || i ==1)
        omega(j) = wroot(i);
        check_omega = omega(j);
        omegaval(j) = fval(i);
        lambda(j) = 2*b/(1+b^2*omega(j)^2);
        if(j >1)
            step = omega(j) - omega(j-1);
        else
           step = pi;
        end
        
        j = j+1;
        
    end
    
end
% % figure(1)
% % plot(wroot,fval)
% % figure(2)
% % plot(omega,omegaval)

lambda 
omega 
figure(3)
plot(1:length(lambda),lambda,'*')
ax = gca;
set(gca, 'YScale', 'log','ytick', [0.0001 0.001 0.01 0.1 1]);

% Finding the eigen functions
x = 0:0.01:1;
eigfun = zeros(length(omega),length(x));

figure(6)
for k= 1: length(omega)
    
    if(mod(k,2)~=0)
        % Odd
        for xi = 1:length(x)
          eigfun(k,xi) = (cos(omega(k)*(x(xi)-0.5)))/(sqrt( 0.5 + ( sin(omega(k))/(2*omega(k)) ) ) );
%         eigfun(k,xi) = (sin(omega(k)*(x(xi)-0.5)))/(sqrt(0.5 - ( sin(omega(k))/(2*omega(k)) ) ) );
        end
    else
        %even
        for xi = 1:length(x)
%         eigfun(k,xi) = (cos(omega(k)*(x(xi)-0.5)))/(sqrt( 0.5 + ( sin(omega(k))/(2*omega(k)) ) ) ); 
          eigfun(k,xi) = (sin(omega(k)*(x(xi)-0.5)))/(sqrt(0.5 - ( sin(omega(k))/(2*omega(k)) ) ) );
        end
    end
    
    
%     eigfunval(k,:) = sqrt(lambda(k)).*eigfun(k,:);
    eigfunval(k,:) = eigfun(k,:);
    plot(x,eigfunval(k,:))
    hold on
end







