
% The Stochastic Logistic Model with Gaussian Noise
%dX=f(x)dt + g(x)DB(t)
%f(x)=x(1-x); g(x)=lambda*x
clear all;

n = 1000;   % number of discretization points
T = 100;    % length of simulation interval
dt = T/n;   % size of time step
randn('state',0); % sets the seed of the random number generator
% simulate a trajectory
x_euler = zeros(n+1,1); % initialization of the trajectory
x0 = 1;        % initial condition
x_euler(1) = x0;        % the initial condition
 
sigma = 0.1;  % parameter values
k1=0.8; r=0.5; k2=1.2; eps=1; S= 0.5; theta=1; gamma=5; phi=(k1+k2)/2;
k=@(x)k1+(k2-k1)/(1+exp(-gamma*(x-phi)));

for j=2:n+1             % the Euler-Maruyama scheme
    dW = sqrt(dt)*randn; % the Wiener increment
    x = x_euler(j-1);
    %x_euler(j) = x + ((x.*(b-d*exp(-a*x))))*dt + sigma*x*dW;  % Ref. D5
   %x_euler(j) = x + (r*(x.*(1-x./k(x))))*dt +sigma*x*dW;  % Prof. Branna model
    x_euler(j) = x + (-r*(x.*(1-x/S).*(1-x./k(x))))*dt +sigma*x*dW;
end
K=k(x);
mu=@ (x)(r*(x.*(1-x./k(x))));% Deterministic Log
%mu=@ (x)(r*(x.*(x/S-1).*(1-x./k(x))));% Deterministic Allee
Mu=mu(x);
% plot the approximation
plot([0:dt:T],x_euler)
hold on
xlabel('t')
%title('prof. Brannan model with Threshold')
