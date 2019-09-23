clear all; close all;clc;

x = 0.66164
%%AL
E6 = 7.802E-02;
E8 = 6.841E-02;

a = (E8-E6)/(0.8-0.6);
b = E6-a*0.6;
EAL = (a*x+b)*10^2

mys  = [9.276,8.445,7.802,6.841,6.146];
Es   = [4,5,6,8,10].*0.1;

plot(Es,mys,'.')
hold on 
[beta] = nlinfit(Es,mys,@(beta,x) beta(1).*x.^2+beta(2).*x+beta(3),[1,1,1]);
 plot(Es,beta(1).*Es.^2+beta(2).*Es+beta(3))

% [beta] = nlinfit(Es,mys,@(beta,x) beta(1).*exp(-x.*beta(2))+beta(3),[10,0.1,5]);
% plot(Es,beta(1).*exp(-Es.*beta(2))+beta(3))
EAL = beta(1).*x.^2+beta(2).*x+beta(3)
%% Lead
E6 = 1.248E-01;
E8 = 8.870E-02;

a = (E8-E6)/(0.8-0.6);
b = E6-a*0.6;
Epb = (a*x+b)*10^2

%% CU
E6 = 7.625E-02;
E8 = 6.605E-02;

a = (E8-E6)/(0.8-0.6);
b = E6-a*0.6;
Ecu = (a*x+b)*10^2

%% Brass
