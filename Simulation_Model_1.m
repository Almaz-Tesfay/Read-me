
% The Stochastic Logistic Model with Gaussian Noise
%dX=f(x)dt + g(x)DB(t)
%f(x)=x(1-x); g(x)=lambda*x
clear all;

n = 1000;   % number of discretization points
T = 20;    % length of simulation interval
dt = T/n;   % size of time step
randn('state',0); % sets the seed of the random number generator
% simulate a trajectory
x_euler = zeros(n+1,1); % initialization of the trajectory
x0 = 1;        % initial condition
x_euler(1) = x0;        % the initial condition
 
sigma = 0.5;  % parameter values

p=1; q=2; mu=q/p; s=0.8;% q>p>0 and s>0. Parameters
f=@(x) -p*x.*(1-mu*exp(-s*x));
for j=2:n+1             % the Euler-Maruyama scheme
    dW = sqrt(dt)*randn; % the Wiener increment
    x = x_euler(j-1);
    %x_euler(j) = x + ((x.*(b-d*exp(-a*x))))*dt + sigma*x*dW;  % Ref. D5
   %x_euler(j) = x + (r*(x.*(1-x./k(x))))*dt +sigma*x*dW;  % Prof. Branna model
    x_euler(j) = x -p*x.*(1-mu*exp(-s*x))*dt +sigma*x*dW;
end

mu=@ (x)-p*x.*(1-mu*exp(-s*x));% Deterministic Log
%mu=@ (x)(r*(x.*(x/S-1).*(1-x./k(x))));% Deterministic Allee
Mu=mu(x);
% plot the approximation
plot([0:dt:T],x_euler)
hold on
xlabel('t')
%title('prof. Brannan model with Threshold')
