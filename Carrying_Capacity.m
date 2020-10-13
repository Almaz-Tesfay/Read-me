%dX = rX(1-X/T)(1-X/K(x))-lambdaX.,  beta=K/T,
% State dependet carrying capacity
clear all;
x1=0; x2=2; dx=1/2^8; xx=x1:dx:x2;
r=0.5; S=0.2;   % parameter values
k1=0.8; k2=1.2;    phi=(k1+k2)/2; gamma=15;% gamma must be near gamC
V=@(x)k1+(k2-k1)/(1+exp(-gamma*(x-phi)));  gamC=4/(k2-k1); % bifurcation point 

M=round((x2-x1)/dx);
X=zeros(1,M+1);
for i=1:M+1
    X(i)=V(xx(i));   % For Carrying capacity
end
plot(xx,X,'.')
hold on
xlabel('X')
ylabel('M(X,\beta)')
