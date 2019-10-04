clear all; close all;clc;

x = 0.66164
%%AL
E6 = 7.802E-02;
E8 = 6.841E-02;

a = (E8-E6)/(0.8-0.6);
b = E6-a*0.6;
EAL = (a*x+b)*10^2

mys  = [9.276,8.445,7.802,6.841,6.146].*1e-02;
Es   = [4,5,6,8,10].*0.1;

mys  = [8.445,7.802,6.841,6.146].*1e-02;
Es   = [5,6,8,10].*0.1;


% plot(Es,mys,'.')
% hold on 
% [beta] = nlinfit(Es,mys,@(beta,x) beta(1).*x.^2+beta(2).*x+beta(3),[1,1,1]);
%  plot(Es,beta(1).*Es.^2+beta(2).*Es+beta(3))

% [beta] = nlinfit(Es,mys,@(beta,x) beta(1).*exp(-x.*beta(2))+beta(3),[10,0.1,5]);
% plot(Es,beta(1).*exp(-Es.*beta(2))+beta(3))

% EAL = beta(1).*x.^2+beta(2).*x+beta(3)


fitfun = @(beta,x) beta(1)*exp(beta(2)*x) + beta(3);
% fitfun = @(beta,x) beta(1).*x.^2+beta(2).*x+beta(3)

[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(Es,mys,@(beta,x)fitfun(beta,x),[1,1,1]);
figure
xs = linspace(0.3,1.3,100);
[Ypred,delta] = nlpredci(@(beta,x)fitfun(beta,x),xs,beta,R,'jacobian',J,'alpha',0.35);
plot(xs,Ypred)
hold on 
plot(Es,mys,'.')

[Ypred,delta] = nlpredci(@(beta,x)fitfun(beta,x),x,beta,R,'jacobian',J,'alpha',0.35);
EAL = 2.699e+00*Ypred
sigmaAL = 2.699e+00*delta
gammelPredAL = (a*x+b)*10^2*2.699e+00



%% Lead
E6 = 1.248E-01;
E8 = 8.870E-02;

a = (E8-E6)/(0.8-0.6);
b = E6-a*0.6;
Epb = (a*x+b)*10^2



mys  = [2.323E-01,1.614E-01,1.248E-01,8.870E-02,7.102E-02,5.876E-02];
Es   = [4,5,6,8,10,12.5].*0.1;
mys  = [1.614E-01,1.248E-01,8.870E-02,7.102E-02];
Es   = [5,6,8,10].*0.1;

fitfun = @(beta,x) beta(1)*exp(beta(2)*x) + beta(3);
% fitfun = @(beta,x) beta(1).*x.^2+beta(2).*x+beta(3)

[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(Es,mys,@(beta,x)fitfun(beta,x),[1,1,1]);
figure
xs = linspace(0.3,1.3,100);
[Ypred,delta] = nlpredci(@(beta,x)fitfun(beta,x),xs,beta,R,'jacobian',J,'alpha',0.35);
plot(xs,Ypred)
hold on 
plot(Es,mys,'.')

[Ypred,delta] = nlpredci(@(beta,x)fitfun(beta,x),x,beta,R,'jacobian',J,'alpha',0.35);
Epb = 1.135E+01*Ypred*10^2
sigmaEpb = 1.135E+01*delta*10^2
gammelPred = (a*x+b)*10^2*1.135E+01




%% CU
E6 = 7.625E-02;
E8 = 6.605E-02;




a = (E8-E6)/(0.8-0.6);
b = E6-a*0.6;
Ecu = (a*x+b)*10^2


mys  = [8.362E-02,7.625E-02,6.605E-02,5.901E-02];
Es   = [5,6,8,10].*0.1;

fitfun = @(beta,x) beta(1)*exp(beta(2)*x) + beta(3);
% fitfun = @(beta,x) beta(1).*x.^2+beta(2).*x+beta(3)

[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(Es,mys,@(beta,x)fitfun(beta,x),[1,1,1]);
figure
xs = linspace(0.3,1.3,100);
[Ypred,delta] = nlpredci(@(beta,x)fitfun(beta,x),xs,beta,R,'jacobian',J,'alpha',0.35);
plot(xs,Ypred)
hold on 
plot(Es,mys,'.')

[Ypred,delta] = nlpredci(@(beta,x)fitfun(beta,x),x,beta,R,'jacobian',J,'alpha',0.35);
ECu = 8.960e+00*Ypred*10^2
sigmaECu = 8.960e+00*delta*10^2
gammelPredCu = (a*x+b)*10^2*8.960e+00



%% Brass
