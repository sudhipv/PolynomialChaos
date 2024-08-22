%% Code to get lambdas, omegas and multipliers for exponential covariance
% For a square domain 
% 0 1

dbounds = zeros(2,2);

prompt = 'Provide domain length along 1st dimension: Eg: 1 for unit-square ';
dbounds(1,2) = input(prompt);

% prompt = 'Please provide mesh length along 2nd dimension: ';
% dbounds(2,2) = input(prompt);

xlims = dbounds(1,:);
%ylims = dbounds(2,:);

a = (xlims(2)-xlims(1))/2;                
b = 1.0; 
x = linspace(xlims(1),xlims(2),100);
x_centered = x-xlims(1)-a;
offsets(1) = xlims(1)+a;

wn = 0;
omegas = [];
lambdas = [];
multipliers = [];

for n = 1:200

%%%%%%%%% Odd terms in x
f = @(w)(1/b - w*tan(a*w));        %% Function to find w

wn = wn+1;
w = fzero(f, wn);                  %% Root of nonlinear function
lambda = 2*b/(1+b^2*w^2);          %% Eigen values
omegas = [omegas ; w];
lambdas = [lambdas ; lambda];
disp(['1st term in x: w = ', num2str(w), ' , lambda = ', num2str(lambda)]);

end 
plot(w,fval)
for n = 1:200

%%%%%%%%% Even terms in x
f = @(w)(w + (1/b)*tan(a*w));

wn = wn+1;
w = fzero(f, wn);   
if w > 1
lambda = 2*b/(1+b^2*w^2);
omegas = [omegas ; w];
lambdas = [lambdas ; lambda];
disp(['2nd term in x: w = ', num2str(w), ' , lambda = ', num2str(lambda)]);
end 

end 
