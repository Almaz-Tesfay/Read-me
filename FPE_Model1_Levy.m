% compute pdf of fpe corr to Levy noise with absorbing BC 
% in a bounded, symmetric domain (-r,r);
% central differencing for 2nd-order deriv 
% and using one-sided differencing near the boundaries;
% 3rd-order WENO for 1st-order deriv
% 3rd-order TVD RK in time
% dx(t)=g*f(x)dB(t)+sigma*Sigma(x)dL(t)

%%% Dynamical analysis on twos to chastic single-species models 
% function [X,U]=tgao7_Dwell(Alp,T)

tic;
clear
alpha=1.5;    T=50;          %%%%change
  
 eps=0.2;
% cigma=@(x)  abs(x).^alpha.*sign(x);
eps=eps^alpha;

d=0;
x0=0.5;
s=1.6; p=0.5; q=1;  mu=2;  x1=1/s*(log(mu)); %mu=q/p; ;% gb>gr>0 and ga>0. Parameters
f=@(x) -p*x.*(1-mu*exp(-s*x)); %  drift term
cigma=@(x)  eps*abs(x).^alpha;% noise term
rr=1;             %right
rl=1;             %left
fpmax=abs(f(-rr));

h=1/45;
J=rr/h;             %right
L=rl/h;             %left
dt = 0.05*h^alpha; 
dtdx = dt/h;
%T=2;

Jt=L+J-1;     %Jt=total number of unknowns  
x=-(rl+rr):h:(rl+rr);

cons=alpha*gamma((1+alpha)/2)/(2^(1-alpha)*sqrt(pi)*gamma(1-alpha/2));
C=-zeta(alpha-1)*h^(-alpha)*cons*eps;          % correction term u''(x)
% coeff for the diffusion and the correction term
Chh=( d/2  )/h^2;
c1=eps*cons/alpha;
c2=eps*cons*h;

b=zeros(Jt,1);
a=zeros(Jt,1);
c=zeros(Jt,1); 

%nonintegral part

% coefficient of U_j
b(2:Jt-1) = (-2*Chh -2*C*cigma(x(J+3:2*J+L-1))- c1*cigma(x(J+3:2*J+L-1)).*(1./(x(J+3:2*J+L-1)+rl).^alpha+1./(rr-x(J+3:2*J+L-1)).^alpha))';
% one-sided diff near boundaries
b(1)  = 2*Chh+2*C*cigma(x(J+2))- c1*cigma(x(J+2)).*(1/(rl+x(J+2))^alpha+1/(rr-x(J+2))^alpha); 
b(Jt) = 2*Chh+2*C*cigma(x(2*J+L)) - c1*cigma(x(2*J+L)).*(1/(rl+x(2*J+L))^alpha+1/(rr-x(2*J+L))^alpha); % one-sided diff
a= Chh*ones(Jt,1)+C*cigma(x(J+1:2*J+L-1)).';  % coefficient of U_(j-1)
c= Chh*ones(Jt,1)+C*cigma(x(J+3:2*J+L+1)).';  % coefficient of U_(j+1) 
c(1)  = -5*Chh-5*C*cigma(x(J+3)); % one-sided diff
a(Jt) = -5*Chh-5*C*cigma(x(2*J+L-1)); % one-sided diff
vp2 = zeros(Jt,1); vp2(3) = 4*Chh+4*C*cigma(x(J+4));  % one-sided diff
vp3 = zeros(Jt,1); vp3(4) =  -Chh-C*cigma(x(J+5)); % one-sided diff 
vm2 = zeros(Jt,1); vm2(Jt-2) = 4*Chh+4*C*cigma(x(2*J+L-2));  % one-sided diff
vm3 = zeros(Jt,1); vm3(Jt-3) =  -Chh-C*cigma(x(2*J+L-3)); % one-sided diff 

% integral part
for j=-L+1:J-1
   b(j+L)= b(j+L) - c2*cigma(x(J+L+1+j))*( sum(1./abs(x(J+2-j:L+J)).^(1+alpha)) ...
                       + sum(1./abs(x(L+J+2:2*J+L-j)).^(1+alpha)) ...
         + .5/abs(x(J+1-j))^(1+alpha) + .5/abs(x(2*J+L+1-j))^(1+alpha) );  
end 
A=spdiags([vm3 vm2 [a(2:end); 0] b ...
           [0; c(1:end-1)] vp2 vp3],-3:3,Jt,Jt);

% coefficient of u_(j+k) 
B=zeros(size(A));
for j=-L+1:J-1
  B(L+j,:)=cigma(x(J+2:2*J+L)).*[1./abs(x(J+2-j:L+J)).^(1+alpha)  0  1./abs(x(L+J+2:2*J+L-j)).^(1+alpha)];
end


%A = diag(ones(Jt,1))+ dt*(A+c2*B); % the iterative matrix for each time-step

X=-rl+rr:h:rr+rr;
%UU=sqrt(pi)*(r^2-X.^2).^(alpha/2)/(2^alpha*gamma(1+alpha/2)*gamma(1/2+alpha/2)); 

UU=sqrt(40/pi)*exp(-40*(X-x0).^2); %gaussian
%UU = zeros(2*J+1,1)';
%lo = (x0+1)/h+1;
%UU(lo) = 1/(2*h);    UU(lo+1) = 1/(2*h);  UU(lo-1) = 1/(2*h);
%UU=ones(size(X))./(rr+rl);
U=UU(2:end-1)';
Un=U; U1=U; U2=U;
%U=UU(2:end-1);


nft=round(T/dt);

nu = length(U);
data = zeros(nu+4,1);
x_star=zeros(nft+1,1);
x_star(1)=x0;
for nt=1:nft

    U1 = U + dt*(A+c2*B)*U;
    % global Lax-Friedrichs(LF) flux splitting
    data(3:nu+2) = (f(X(2:end-1)) + fpmax)'.*U/2;
    fx1 = derWENOr2_minus(data,h);
    data(3:nu+2) = (f(X(2:end-1)) - fpmax)'.*U/2;
    fx2 = derWENOr2_plus(data,h);
    U1 = U1 - dtdx*(fx1+fx2);
    
    U2 = 0.75*U + U1/4 + (dt/4)*(A+c2*B)*U1;
    data(3:nu+2) = (f(X(2:end-1)) + fpmax)'.*U1/2;
    fx1 = derWENOr2_minus(data,h);
    data(3:nu+2) = (f(X(2:end-1)) - fpmax)'.*U1/2;
    fx2 = derWENOr2_plus(data,h);
    U2 = U2 - (dtdx/4)*(fx1+fx2);
    
    Un = U/3 + 2*U2/3 + (2*dt/3)*(A+c2*B)*U2;
    data(3:nu+2) = (f(X(2:end-1)) + fpmax)'.*U2/2;
    fx1 = derWENOr2_minus(data,h);
    data(3:nu+2) = (f(X(2:end-1)) - fpmax)'.*U2/2;
    fx2 = derWENOr2_plus(data,h);
    Un = Un - (2*dtdx/3)*(fx1+fx2);     
     U=Un;
  [ind3]=find(U==(max(U)));
    x_star(nt+1)=X(min(ind3)+1);
end

t=0:dt:nft*dt;


 %figure
 %plot(t,x_star)
 %figure
 plot(X(2:end-1)',U,'-')
 xlabel('X')
 %title('T=50,  \alpha=1, r=0.1')
%legend('\sigma=0.1','\sigma=0.4','\sigma=0.7','\sigma=1')
 %legend('\epsilon=0.0', '\epsilon=0.5', '\epsilon=1.0')
 ylabel('pdf-P')
 hold on
 %grid on
 save  A2 alpha  x_star
toc;
