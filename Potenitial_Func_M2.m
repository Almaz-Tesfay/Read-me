%dX = rX(1-X/T)(1-X/K)-lambdaX.,  beta=K/T,
%%0 < lambda/r < ((1+beta)^2/4beta  -1)
clear all;
x1=0; x2=1.5; dx=1/2^6;
xx=x1:dx:x2;
r=0.1; S=0.4;   % parameter values
k1=0.8; r=1; k2=1.2;    phi=(k1+k2)/2; gamma=12;% gamma must be near gamC
V=@(x)k1+(k2-k1)/(1+exp(-gamma*(x-phi)));
gamC=4/(k2-k1) % bifurcation poin 
mu=@ (x)r*(x.*(1-x./V(x)));% Deterministic Logistic
%mu=@ (x)(r*(x.*(x/S-1).*(1-x./V(x))));% Deterministic Allee
%Mu=integral (mu(x))
Mu=@ (x)r*(1-1/4*gamma*(k2-k1)*(x).^2/2);
M=round((x2-x1)/dx);
X=zeros(1,M+1);
for i=1:M+1
 X(i)=Mu(xx(i));
end
plot(xx,X)
hold on
%xlabel('y')
title('Potential Function')
%title('Potential Fun for Logistic Growth with Threshold')
%ylabel('V(y)')
grid on