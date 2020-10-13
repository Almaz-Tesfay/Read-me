%%% Dynamical analysis on twos to chastic single-species models 
tic;
clear all;
x1=0; x2=1; dx=1/2^10;
xx=x1:dx:x2;
s=1;
p=1; q=2; mu= 2; % mu=q/p; % gb>gr>0 and ga>0. Parameters
f=@(x)  -p*x.*(1-mu*exp(-s*x)); %  Dterministic model

M=round((x2-x1)/dx);
X=zeros(1,M+1);
for i=1:M+1
    X(i)=f(xx(i));
end
plot(xx,X,'-')
hold on
%xlabel('y')
%title('Phase diagram')
%title('Potential Fun for Logistic Growth with Threshold')
%ylabel('V(y)')
%grid on
