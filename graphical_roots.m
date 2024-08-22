


%%% Graphical method of finding roots


% Trying to plot the function to find the approximate roots

clear all
w = 4.8:.001:4.9;
b = 2;
a = 1;
val = zeros(1,length(w));


for i = 1:length(w)
    % odd
%     val(i) = ((1/b)-w(i)*tan(w(i)*a));
  
    
    % even
      val(i) = (w(i)+(1/b)*tan(w(i)*a)) ;

      
end

figure(5)
plot(w,zeros(1,length(w)),'-')
hold on
plot(w,val,'b')
hold on 
line(w,0)

% Roots from the graphical method 
%..............................%
% Odd       Even
% .......................
% 1.3065    3.6733
% 6.5846    9.6318
% 12.7234   15.8345
% 18.9550   22.0819
% 25.2120   28.3455
 
%...............................%