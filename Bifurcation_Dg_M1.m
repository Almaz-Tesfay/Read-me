
%%%% Plotting the bifurcation diagram 
%Bifurcation diagram of Model 1
%plots for mu<1
ezplot('x^2-x*log(mu)',[0 1 -1 1])  % mu=q/p; s=1
hold on
%plots for mu>1
ezplot('x^2-x*log(mu)',[1 2 -1 1])  % mu=q/p;s=1

