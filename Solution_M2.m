%dX = rX(1-X/T)(1-X/K(x))-lambdaX.,  beta=K/T,


t=[0:0.01:50];               % vector of time points
 r=0.1;   % parameter values
 k1=0.8; k2=1.2;    phi=(k1+k2)/2; gamma=10;% gamma must be near gamC
% V=@(x)k1+(k2-k1)/(1+exp(-gamma*(x-phi)));
 gamC=4/(k2-k1); % bifurcation poin 

A=-r/gamC.*(gamC-gamma);
 u0=0.2;
 u=u0*exp(A*t); % u(t) is the solution                 % vector of positions
axis([0 1 0 50]);
plot(t,u,'r-.'), hold on         % plot this trajectory, stay on this plot
u0=0.5;
u=u0*exp(A*t);
plot(t,u,'b--'), hold on
u0=0.7;
u=u0*exp(A*t);
plot(t,u,'gr.'), hold on
u0=0.9;
 u=u0*exp(A*t);
plot(t,u,'y-'), hold on
u0=1;
u=u0*exp(A*t);
plot(t,u,'r.-'), hold on
%%grid on
title('Trajectories');
xlabel('time');



% 
% clear all;
% t1=0; t2=5; dt=0.05;
% t=t1:dt:t2;
% r=0.5; S=0.2;   % parameter values
% k1=0.8; k2=1.2;    phi=(k1+k2)/2; gamma=12;% gamma must be near gamC
% V=@(x)k1+(k2-k1)/(1+exp(-gamma*(x-phi)));
% gamC=4/(k2-k1) % bifurcation poin 
% % Deterministic Logistic growth model
% A=-r*((1-1/gamC)*gamma)
% u0=0.5;
% u= @ (x) u0*exp(A*t); % u(t) is the solution
% M=round((t2-t1)/dt);
% T=zeros(1,M+1);
% 
% for i=1:M+1
%     T(i)=u(t(i));   % For Carrying capacity
% end
% plot(t,T,'-')
% hold on
% xlabel('t')
% ylabel('u')
