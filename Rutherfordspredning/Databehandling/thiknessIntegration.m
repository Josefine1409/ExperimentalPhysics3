clear all; close all;clc;

c2E = @(x) (x*0.76535+15.78)./1000
Ein = 0.349

mG=196.966-4.4858e-4*79-0.03343120468;
mC=12-4.4858e-4*12;

I = importfile('..\Data\Thickness\Sample2_0degree_480sec.asc');
E = 1:length(I);
Ierr = sqrt(I);
errorbar(E,I,Ierr,'.')

peakBorders = [242,353;353,437];
name = 'Sample1_0degree'
data = fitGaussInSpectrum(E,I,Ierr,name,2,peakBorders);
EOutC = c2E(data(1,1))
EOutCErr = c2E(data(2,1))
EOutAu = c2E(data(1,2))
EOutAuErr = c2E(data(2,2))

xs = linspace(1e-8,3e-7,20);

% Esfes =linspace(0.0,1,100);
% figure
% plot(Esfes,stoppingpowerC(Esfes))
% Eouts= 0;
% for i = 1:length(xs)
%     Eouts(i) = EStopTop(Ein,xs(i),pi-160/180*pi,@(E) stoppingpowerC(E),mC);
% end
% figure
% plot(xs,Eouts)
% xC = fzero(@(x) (EStopTop(Ein,x,pi-160/180*pi,@(E) stoppingpowerC(E),mC)-EOutC) , 1e-7)


xs = linspace(1e-10,1e-8,20);

% Esfes =linspace(0.0,1,100);
% figure
% plot(Esfes,stoppingpowerC(Esfes))
Eouts= 0;
for i = 1:length(xs)
    Eouts(i) = EStopBottom(Ein,2.5e-8,xs(i),pi-160/180*pi,@(E) stoppingpowerC(E),@(E) stoppingpowerAu(E),mG);
end
figure
plot(xs,Eouts)



% xAu = fzero(@(x) (EStopBottom(Ein,1.5e-7,x,pi-160/180*pi,@(E) stoppingpowerC(E),@(E) stoppingpowerAu(E),mG)-EOutAu) , 1e-7)





function data = fitGaussInSpectrum(X,Y,Yerr,name,n,peakBorder)


figure
histogram('BinEdges',[0,X(1:end)],'BinCounts',Y,'EdgeColor','none')
hold on
xlabel('Channel number (Ch)')
ylabel('Counts (n)')
set(gca,'FontSize',15) 
title(name)

for i = 1:n

x2 =   peakBorder(i,1);
x3 =   peakBorder(i,2);  

x = X;
y = Y;
yerr = Yerr;
higherIndex = x>x2;
x = x(higherIndex);
y = y(higherIndex);
yerr = yerr(higherIndex);

lowerIndex = x<x3;
x = x(lowerIndex);
y = y(lowerIndex);
yerr = yerr(lowerIndex);

beta0 = [(x2+x3)/2,(x3-x2)/3,max(y)/2,0,20];

w = 1./(yerr+1).^2;
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x',y,@fitfunction,beta0,'weights',w);
plot(x,fitfunction(beta,x),'k--','linewidth',2)
us = CovB/MSE;
mse =MSE;
MSECount(i) = MSE;
pValue(i) = 1-chi2cdf(MSE*(length(y)-5),(length(y)-5));
P_Value = 1-chi2cdf(MSE*(length(y)-5),(length(y)-5));

peakChannel(i) = beta(1,1);
peakUns(i) = us(1,1);
ci = nlparci(beta,R,'jacobian',J,'alpha',0.35);
peakUns(i) = (ci(1,2)-ci(1,1))/2;

end

peakUns;
peakChannel;
data = [peakChannel;peakUns;MSECount];
end


function y = fitfunction(beta,x)
    y = beta(4).*x+beta(5)+beta(3).*exp(-((x-beta(1))./(beta(2))).^2./2);
end




function Eout=EStopTop(Ein,x,theta,S,m)
m1=1.007276;

[q E]=ode45(@(x,E) -S(E),[0 x/2],Ein);
E1=E(end);
K2=@(theta,m2)((m1*cos(theta)+sqrt(m2^2-m1^2*sin(theta).^2))./(m1+m2)).^2;
[q E]=ode45(@(x,E) -S(E),[0 x/(2*cos(theta))],E1*K2(pi-theta,m));

Eout=E(end);
end

function Eout=EStopBottom(Ein,x1,x2,theta,S1,S2,m)
mH=1.007276;

[x E]=ode45(@(x,E) -S1(E),[0 x1],Ein);
E1=E(end);
[x E]=ode45(@(x,E) -S2(E),[0 x2/2],E1);
E1=E(end);
K2=@(theta,m2)((mH*cos(theta)+sqrt(m2^2-mH^2*sin(theta).^2))./(mH+m2)).^2;
[x E]=ode45(@(x,E) -S2(E),[0 x2/(2*cos(theta))],E1*K2(pi-theta,m));
E1=E(end);
[x E]=ode45(@(x,E) -S1(E),[0 x1/(cos(theta))],E1);

Eout=E(end);
end