clear all;
close all;
x1=0; x2=2; dx=1/2^6;
xx=x1:dx:x2;
% Parameters q>p, s >0

mu = [0.5 1 2];

s=1; p=1; 

figure
hold on

for K = 1:3
q = mu(K)*p;
V2=@(x)exp(-s*x);
V1=@ (x)(p/2)*x.^2+q*V2(x).*(s*x+1)/s^2;
M=round((x2-x1)/dx);
X=zeros(1,M+1);
for i=1:M+1
 X(i)=V1(xx(i));
end
plot(xx,X,'-','LineWidth',2)

end
h2 = legend(strcat('$\mu=$',num2str(mu(1))),strcat('$\mu=$',num2str(mu(2))),strcat('$\mu=$',num2str(mu(3))),'Location','best');
set(h2,'Interpreter','latex')
plot([0,0],ylim,'k-')
xlabel('x')
title('Potential Function')
box on
