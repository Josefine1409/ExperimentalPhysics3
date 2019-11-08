% clear all; close all;clc;
    
Ein = 349
Ein = 349.9

c2E = @(x) x.*0.76675+11.5073

m1=1.007276;
mG=196.966-4.4858e-4*79-0.03343120468;
mC=12-4.4858e-4*12;

K2=@(theta,m2)((m1*cos(theta)+sqrt(m2^2-m1^2*sin(theta).^2))./(m1+m2)).^2;


theta =       [30,            40];
peakValues =  {c2E([434,445]),c2E([434,431])};
peakBorders = {[390,480],     [380,490]};


thetas=theta./180.*pi;

for i = 1:length(thetas)
    peakValues{i} = [EoutNy(Ein,thetas(i),1),EoutNy(Ein,thetas(i),2)];
end


linFun =@(beta,x) (x-beta(2))/beta(1);
dataAg = [];
dataC = [];

addpath('\ipf13')

Res ={}

i = 2

[X,Y,Yerr] = hisFraData(['..\Data\AngularDependency\' num2str(theta(i)) 'degree.asc']);
W = 1/Yerr.^2;
peakValue = peakValues{i}
peakChannel = (peakValues{i}-11.5073)./0.76675;
pBorders = peakBorders{i}
Xs = pBorders(1):pBorders(2)
Ys = Y(Xs)

pb = peakBorders{i}
start = [peakChannel(1),14.4,peakChannel(2),16.8]

figure
histogram('BinEdges',[1/2,Xs(1:end)+1/2],'BinCounts',Ys,'EdgeColor','none')

hold on 

x01 =       437.5 % (437.2, 437.7)
x02 =       434.1  %(433.8, 434.4)
a1 =   1.092e+04  %(8833, 1.3e+04)
a2 =   1.726e+04  %(1.52e+04, 1.933e+04)
sigma1 =       4.491 % (4.271, 4.711)
sigma2 =       6.359 % (6.249, 6.469)
c1 =      -1.733  %(-2.366, -1.1)
c2 =       948.4  %(671.2, 1226)

EAu2 = c2E(x01)
EC2 = c2E(x02)


E0Au2Err = -(437.2-437.7)/2*0.76675
E0C2Err = -(433.8-434.4)/2*0.76675



a1Err =  (1.3e+04-8833)/2
sigma1Err = (4.711-4.271)/2

counts2Au = 2*sqrt(pi)*a1*sigma1

countUns2Ai = sqrt(...
    (2*sqrt(pi)*a1).^2*sigma1Err.^2+...
    (2*sqrt(pi)*sigma1).^2*a1Err^2)

a2Err =  (1.933e+04-(1.52e+04))/2
sigma2Err = (6.469-6.249)/2

counts2C = 2*sqrt(pi)*a2*sigma2

countUns2C = sqrt(...
    (2*sqrt(pi)*a2).^2*sigma2Err.^2+...
    (2*sqrt(pi)*sigma2).^2*a2Err^2)


plot(Xs,c1*Xs+c2+a1.*exp(-((Xs-x01)./(2*sigma1)).^2)+a2.*exp(-((Xs-x02)./(2*sigma2)).^2),'-.','linewidth',3)
plot(Xs,c1*Xs+c2+a1.*exp(-((Xs-x01)./(2*sigma1)).^2),'--','linewidth',1)
plot(Xs,c1*Xs+c2+a2.*exp(-((Xs-x02)./(2*sigma2)).^2),'--','linewidth',1)
plot(Xs,c1*Xs+c2,'--','linewidth',1)
legend('Fit','Data','Fit confidence','Location','southwest')

plot([peakChannel(1),peakChannel(1)],[0,max(Ys)],'k')
plot([peakChannel(2),peakChannel(2)],[0,max(Ys)],'k')



i = 1

[X,Y,Yerr] = hisFraData(['..\Data\AngularDependency\' num2str(theta(i)) 'degree.asc']);
W = 1/Yerr.^2;
peakValue = peakValues{i};
peakChannel = (peakValues{i}-11.5073)./0.76675;
pBorders = peakBorders{i};
Xs = pBorders(1):pBorders(2);
Ys = Y(Xs);

pb = peakBorders{i};
start = [peakChannel(1),14.4,peakChannel(2),16.8];

figure
histogram('BinEdges',[Xs(1)-1/2,Xs(1:end)+1/2],'BinCounts',Ys,'EdgeColor','none')

hold on 

% x01 =       434.7  %(434.7, 434.7)
% x02 =       439.6  %(439, 440.2)
% a1 =   8.655e+04  %(8.539e+04, 8.771e+04)
% a2 =   1.236e+04  %(1.128e+04, 1.343e+04)
% sigma1 =       5.069  %(7.118, 7.22)
% sigma2 =       10.45  %(14.09, 15.45)
% c1 =      -3.389  %(-6.988, 0.2111)
% c2 =        1742  %(268.4, 3216)


% x01 =       436.2  %(435.8, 436.5)
% x02 =       434.4  %(434.2, 434.5)
% a1 =   3.844e+04  %(3.426e+04, 4.262e+04)
% a2 =   5.959e+04  %(5.557e+04, 6.362e+04)
% sigma1 =           7%  (fixed at bound)
% sigma2 =       4.585 % (4.442, 4.728)
% c1 =       12.94  %(6.945, 18.94)
% c2 =       -4387  %(-7008, -1765)

x01 =       435.4 % (435, 435.7)
x02 =         434  %(433.4, 434.6)
a1 =   7.225e+04  %(4.944e+04, 9.506e+04)
a2 =   2.526e+04  %(2594, 4.793e+04)
sigma1 =           6%  (5.597, 6.403)
sigma2 =       3.909 % (2.976, 4.842)
c1 =       17.94  %(9.353, 26.53)
c2 =       -6169  %(-9900, -2438)

EAu1 = c2E(x01)
EC1 = c2E(x02)
E0Au1Err = -(435- 435.7)/2*0.76675
E0C1Err = -(433.4-434.6)/2*0.76675


a1Err =  -(3.426e+04-4.262e+04)/2
sigma1Err = -(5.597- 6.403)/2

counts1Au = 2*sqrt(pi)*a1*sigma1

countUns1Ai = sqrt(...
    (2*sqrt(pi)*a1).^2*sigma1Err.^2+...
    (2*sqrt(pi)*sigma1).^2*a1Err^2)

a2Err =  -(2594- 4.793e+04)/2
sigma2Err = -(2.976-4.842)/2

counts1C = 2*sqrt(pi)*a2*sigma2

countUns1C = sqrt(...
    (2*sqrt(pi)*a2).^2*sigma2Err.^2+...
    (2*sqrt(pi)*sigma2).^2*a2Err^2)






plot(Xs,c1*Xs+c2+a1.*exp(-((Xs-x01)./(2*sigma1)).^2)+a2.*exp(-((Xs-x02)./(2*sigma2)).^2),'-.','linewidth',3)
plot(Xs,c1*Xs+c2+a1.*exp(-((Xs-x01)./(2*sigma1)).^2),'--','linewidth',1)
plot(Xs,c1*Xs+c2+a2.*exp(-((Xs-x02)./(2*sigma2)).^2),'--','linewidth',1)
plot(Xs,c1*Xs+c2,'--','linewidth',1)
legend('Fit','Data','Fit confidence','Location','southwest')

plot([peakChannel(1),peakChannel(1)],[0,max(Ys)],'k')
plot([peakChannel(2),peakChannel(2)],[0,max(Ys)],'k')
















function [X,Y,Yerr] = hisFraData(filename)
filename;
spectrum = importfile(filename);
X = 1:length(spectrum);
Y = spectrum;
Yerr = sqrt(Y);
end
function eout = EoutNy(Ein,Theta,id)
Ein = Ein/1000
mG=196.966-4.4858e-4*79-0.03343120468;
mC=12-4.4858e-4*12;
m1=1.007276;

TAu = 25e-10;
TC = 250e-10;
eout =[];
m = [mG,mC];
K2=@(theta)((m1*cos(theta)+sqrt(m(id)^2-m1^2*sin(theta).^2))./(m1+m(id))).^2;
for theta = Theta;
SC = stoppingpowerC(Ein);
SAu = stoppingpowerAu(Ein);
if (id==1)
    Escatter = (Ein-TAu/2*SAu).*K2(theta);
    SC = stoppingpowerC(Escatter);
    SAu = stoppingpowerAu(Escatter);
    if (theta<pi/2)
        eout(end+1) = Escatter-(TAu./2.*SAu+TC.*SC)./norm(cos(theta));
    else
        eout(end+1) = Escatter-(TAu./2.*SAu)./norm(cos(theta));
    end
else
    Escatter = (Ein-TAu*SAu-TC*SC).*K2(theta);
    SC = stoppingpowerC(Escatter);
    SAu = stoppingpowerAu(Escatter);
    if (theta<pi/2)
        eout(end+1) = Escatter-(TC/2.*SC)./norm(cos(theta));
    else
        eout(end+1) = Escatter-(TC/2.*SC+TAu.*SAu)./norm(cos(theta));
    end
end
end
eout = eout*1000;
end
