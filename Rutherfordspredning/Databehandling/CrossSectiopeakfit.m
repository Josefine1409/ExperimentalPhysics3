clear all; close all;clc;
    
Ein = 349
Ein = 349.9

c2E = @(x) x.*0.76468+12.393
m1=1.007276;
mG=196.966-4.4858e-4*79-0.03343120468;
mC=12-4.4858e-4*12;

K2=@(theta,m2)((m1*cos(theta)+sqrt(m2^2-m1^2*sin(theta).^2))./(m1+m2)).^2;


theta =       [30,            40,           50,             60,            70,            75,            110,           120,           130,           140,           150,           160];
peakValues =  {c2E([434,445]),c2E([434,431]),c2E([433,427]),c2E([433,410]),c2E([433,396]),c2E([428,388]),c2E([432,350]),c2E([433,340]),c2E([433,330]),c2E([433,323]),c2E([433,318]),c2E([433,314])};
peakBorders = {[390,480],     [380,490],     [380,480],     [350,460],     [340,460],     [340,455],     [310,455],     [300,455],     [280,455],     [280,455],     [280,455],     [280,455]};


thetas=[30:10:70 75 110:10:160]./180.*pi;

for i = 1:length(thetas)
    peakValues{i} = [EoutNy(Ein,thetas(i),1),EoutNy(Ein,thetas(i),2)];
end


linFun =@(beta,x) (x-beta(2))/beta(1);
dataAg = [];
dataC = [];

addpath('\ipf13')

Res ={}

for i = 1:length(theta)
    [X,Y,Yerr] = hisFraData(['..\Data\AngularDependency\' num2str(theta(i)) 'degree.asc']);
    peakValue = peakValues{i}
    peakChannel = (peakValues{i}-12.393)./0.76468;

    pb = peakBorders{i}
    start = [peakChannel(1),14.4,peakChannel(2),16.8]
%     (signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,autozero,fixedparameters,plots,bipolar,minwidth,DELTA,clipheight)
    signal = [X;Y'];
    center = (pb(1)+pb(2))/2;
    window= (-pb(1)+pb(2));
    NumPeaks = 2;
    peakshape = 1
    extra = ''
    NumTrials = 5;
    BaselineMode = 1;
    fixedparameters = [];
    plots= 0;
    bipolar= 0;
    minwidth= 12;
    
    [FitResults,GOF,baseline,coeff,residual,xi,yi,bootstrap]= peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials, start, BaselineMode, fixedparameters, plots, bipolar, minwidth);
        
% [FitResults,GOF,baseline,coeff,residual,xi,yi,bootstrap]= peakfit([X;Y'],(pb(1)+pb(2))/2,(-pb(1)+pb(2)),2,1,'',2,start,1,[],1,0,7);
    FitResults
    
    EnergyC(i) = c2E(FitResults(2,2));
    EnergyUnsC(i) = (bootstrap(2,2+4))*0.76468;
    CountsC(i) = (FitResults(2,5));
    CountsUnsC(i) = (bootstrap(2,5+4));
    EnergyAg(i) = c2E(FitResults(1,2));
    EnergyUnsAg(i) = (bootstrap(1,2+4))*0.76468;
    CountsAg(i) = (FitResults(1,5));
    CountsUnsAg(i) = (bootstrap(1,5+4));
    
    
    
figure
hold on
xlabel('Channel number (Ch)')
ylabel('Counts (n)')
set(gca,'FontSize',15) 
% title('Hej')

x2 =   pb(1);
x3 =   pb(2);  
xlim([c2E(x2),c2E(x3)])
% 
x = X;
y = Y;

yerr = Yerr+(Yerr==0);
indexes  = ((x>x2)&(x<x3));
x = x(indexes);
y = y(indexes);
yerr = yerr(indexes);
histogram('BinEdges',c2E([0,X(1:end)]),'BinCounts',Y,'EdgeColor','none')

beta = FitResults;
plot(c2E(x),baseline(1).*x+baseline(2)+beta(1,3).*exp(-((x-beta(1,2))./(beta(1,4))).^2*log(2)*4)+beta(2,3).*exp(-((x-beta(2,2))./(beta(2,4))).^2*log(2)*4),'-.','linewidth',3)
plot(c2E(x),baseline(1).*x+baseline(2)+beta(1,3).*exp(-((x-beta(1,2))./(beta(1,4))).^2*log(2)*4),'--','linewidth',1)
plot(c2E(x),baseline(1).*x+baseline(2)+beta(2,3).*exp(-((x-beta(2,2))./(beta(2,4))).^2*log(2)*4),'--','linewidth',1)
plot(c2E(x),baseline(1).*x+baseline(2),'--','linewidth',1)
legend('Data','Total Fit','Ag peak','C peak','Liniar background')


end

%% Maltekode

ts=[300 300 300 300 350 600 600 600 600 600 600 600];
FCs=[46449 51781 60300 65892 81962 34228 103355 80585 101501 106019 102580 37053];

GCs=CountsAg;
CCs=CountsC;
GEs=EnergyAg;
CEs=EnergyC;

% guld = GCs(1);
% GCs(1) = CCs(1);
% CCs(1) = guld;
% 
% guldE = GEs(1);
% GEs(1) = CEs(1);
% CEs(1) = guldE;

sigmaGCs=CountsUnsAg;
sigmaCCs=CountsUnsC;
sigmaGEs=EnergyUnsAg;
sigmaCEs=EnergyUnsC;

% guldS = sigmaGCs(1);
% sigmaGCs(1) = sigmaCCs(1);
% sigmaCCs(1) = guldS;
% 
% guldES = sigmaGEs(1);
% sigmaGEs(1) = sigmaCEs(1);
% sigmaCEs(1) = guldES;


GCs=GCs./FCs;
CCs=CCs./FCs;
sigmaGCs = sigmaGCs./FCs;
sigmaCCs = sigmaCCs./FCs;


GCs=GCs([3:12])
CCs=CCs([3:1:12])
sigmaGCs = sigmaGCs([3:1:12])
sigmaCCs = sigmaCCs([3:1:12])


cs=@(theta,m,E) 1./(sin(theta/2)).^2.*1./(K2(theta,m).*E).^2;
cs=@(theta,m,E) 1./(sin(theta/2)).^4;
% cs=@(theta,m,E) 1./(sin(theta/2)).^2.*1./(E).^2;


% [betaG,RG,JG,CovB,MSE]=nlinfit(thetas,GCs,@(C,thetas) C(1)*cs(thetas,mG,Ein)+C(2),[1,0]);
% % [betaC,RC,JC,CovB,MSE]=nlinfit(thetas,CCs,@(C,thetas) C(1)*cs(thetas,mC,Ein)+C(2),[1,0]);
% [betaG,RG,JG,CovB,MSE]=nlinfit(thetas([1:5,7:12]),GCs,@(C,thetas) C(1)*cs(thetas,mG,Ein)+C(2),[1,0],'weights',1./sigmaGCs.^2);
% [betaC,RC,JC,CovB,MSE]=nlinfit(thetas([1:5,7:12]),CCs,@(C,thetas) C(1)*cs(thetas,mC,Ein)+C(2),[1,0],'weights',1./sigmaCCs.^2);

[betaG,RG,JG,CovB,MSEG]=nlinfit(thetas([3:1:12]),GCs,@(C,thetas) C(2)+ C(1)*cs(thetas,mG,Ein),[1,0]);
[betaC,RC,JC,CovB,MSEC]=nlinfit(thetas([3:1:12]),CCs,@(C,thetas)C(2)+ C(1)*cs(thetas,mC,Ein),[1,0]);

Chi2TestGold = sum((GCs-betaG(1).*cs(thetas([3:1:12]),mG,Ein)-betaG(2)).^2./(sigmaGCs.^2)) % Chi i anden test
nu = length(GCs)-1 %Antal frihedsgrader
Signifikansniveau = 1-chi2cdf(Chi2TestGold,nu) % P v�rdi
TestN = Chi2TestGold/nu %Ca. test der burde give 1
disp('---------------------')


Chi2TestCarbon = sum((CCs-betaC(1).*cs(thetas([3:1:12]),mC,Ein)-betaC(2)).^2./(sigmaCCs.^2)) % Chi i anden test
nu = length(GCs)-1 %Antal frihedsgrader
Signifikansniveau = 1-chi2cdf(Chi2TestCarbon,nu) % P v�rdi
TestN = Chi2TestCarbon/nu %Ca. test der burde give 1
disp('---------------------')


Thetas=linspace(thetas(3),thetas(end),1000);

fig = figure
errorbar(theta(3:end),GEs(3:end),sigmaGEs(3:end),'r.','markersize',10)
hold on
xlabel('\theta [deg]')
ylabel('Scattering Energy [KeV]')
title('Scattering on Gold')
plot(Thetas/pi*180,EoutNy(Ein,Thetas,1),'linewidth',2)
plot(Thetas/pi*180,Ein*K2(Thetas,mG),'linewidth',2)
set(gca,'FontSize',15) 
xlim([theta(1),theta(end)])

saveas(fig,'EnergyAngleGoldU1p.fig')
saveas(fig,'EnergyAngleGoldU1p.eps','epsc')

fig= figure
errorbar(theta(3:end),CEs(3:end),sigmaCEs(3:end),'r.','markersize',10)
hold on
xlabel('\theta [deg]')
ylabel('Scattering Energy [KeV]')
title('Scattering on Carbon')
plot(Thetas/pi*180,EoutNy(Ein,Thetas,2),'linewidth',2)
plot(Thetas/pi*180,Ein*K2(Thetas,mC),'linewidth',2)
set(gca,'FontSize',15) 
xlim([theta(1),theta(end)])

saveas(fig,'EnergyAngleCarbonU1p.fig')
saveas(fig,'EnergyAngleCarbonU1p.eps','epsc')

fig = figure
[Ypred,deltaY] = nlpredci(@(C,th) C(2)+C(1)*cs(th,mG,Ein),Thetas,betaG,RG,'jacobian',JG,'alpha',0.35);

hold on
xlabel('\theta [Deg]')
ylabel('Normed Counts')
title('Scattering on Gold')
plot(Thetas/pi*180,Ypred,'linewidth',2)
errorbar(theta([3:1:12]),GCs,sigmaGCs,'r.','markersize',10)

plot(Thetas/pi*180,Ypred+deltaY,'k--','linewidth',1)
plot(Thetas/pi*180,Ypred-deltaY,'k--','linewidth',1)
legend('Fit','Data','Fit confidence')
set(gca,'FontSize',15) 
xlim([theta(1),theta(end)])

saveas(fig,'CrossSectionGoldU1p.fig')
saveas(fig,'CrossSectionGoldU1p.eps','epsc')

fig = figure
[Ypred,deltaY] = nlpredci(@(C,th)C(2)+ C(1)*cs(th,mC,Ein),Thetas,betaC,RC,'jacobian',JC,'alpha',0.35);

hold on
xlabel('\theta [deg]')
ylabel('Normed Counts')
title('Scattering on Carbon')
plot(Thetas/pi*180,Ypred,'linewidth',2)
errorbar(theta([3:1:12]),CCs,sigmaCCs,'r.','markersize',10)

plot(Thetas/pi*180,Ypred+deltaY,'k--','linewidth',1)
plot(Thetas/pi*180,Ypred-deltaY,'k--','linewidth',1)
set(gca,'FontSize',15) 
xlim([theta(1),theta(end)])


legend('Fit','Data','Fit confidence')

saveas(fig,'CrossSectionCarbonU1p.fig')
saveas(fig,'CrossSectionCarbonU1p.eps','epsc')

%%

fig = figure
ThetasLinG = 1./(sin(Thetas/2)).^4;
thetasLinG = 1./(sin(theta([3:1:12])/180*pi/2)).^4;

[Ypred,deltaY] = nlpredci(@(C,th)C(2)+ C(1)*th,ThetasLinG,betaG,RG,'jacobian',JG,'alpha',0.35);

hold on
xlabel('1/(sin^4(\theta/2))')
ylabel('Normed Counts')
title('Scattering on Gold')
plot(ThetasLinG,Ypred,'linewidth',2)
errorbar(thetasLinG,GCs,sigmaGCs,'r.','markersize',10)

plot(ThetasLinG,Ypred+deltaY,'k--','linewidth',1)
plot(ThetasLinG,Ypred-deltaY,'k--','linewidth',1)
legend('Fit','Data','Fit confidence','Location','northwest')
set(gca,'FontSize',15) 
xlim([thetasLinG(end),thetasLinG(1)])


saveas(fig,'CrossSectionGoldLinearU1p.fig')
saveas(fig,'CrossSectionGoldLinearU1p.eps','epsc')

fig = figure
ThetasLinC = 1./(sin(Thetas/2)).^4;
thetasLinC = 1./(sin(theta([3:1:12])/180*pi/2)).^4;
[Ypred,deltaY] = nlpredci(@(C,th) C(2)+C(1)*th,ThetasLinC,betaC,RC,'jacobian',JC,'alpha',0.35);

hold on
xlabel('1/(sin^4(\theta/2))')
ylabel('Normed Counts')
title('Scattering on Carbon')
plot(ThetasLinC,Ypred,'linewidth',2)
errorbar(thetasLinC,CCs,sigmaCCs,'r.','markersize',10)

plot(ThetasLinC,Ypred+deltaY,'k--','linewidth',1)
plot(ThetasLinC,Ypred-deltaY,'k--','linewidth',1)
set(gca,'FontSize',15) 
xlim([thetasLinC(end),thetasLinC(1)])


legend('Fit','Data','Fit confidence','Location','northwest')
saveas(fig,'CrossSectionCarbonLinearU1p.fig')
saveas(fig,'CrossSectionCarbonLinearU1p.eps','epsc')


%% Fit
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



function [X,Y,Yerr] = hisFraData(filename)
filename;
spectrum = importfile(filename);
X = 1:length(spectrum);
Y = spectrum;
Yerr = sqrt(Y);
end

%%
% data har sturktur [peakChannel,peakUns,peakValue]

function [data1,data2] = fitGaussInSpectrum(X,Y,Yerr,name,peakValue,peakBorder)
n = length(peakValue);

figure
hold on
xlabel('Channel number (Ch)')
ylabel('Counts (n)')
set(gca,'FontSize',15) 
title(name)

x2 =   peakBorder(1);
x3 =   peakBorder(2);  
xlim([x2,x3])


x = X;
y = Y;

yerr = Yerr+(Yerr==0);

higherIndex = x>x2;
x = x(higherIndex);
y = y(higherIndex);
yerr = yerr(higherIndex);

lowerIndex = x<x3;
x = x(lowerIndex);
y = y(lowerIndex);
yerr = yerr(lowerIndex);

peakChannel = (peakValue-12.393)./0.76468;

beta0 = [0,0,y(round(peakChannel(1))==x),peakChannel(1),5,y(round(peakChannel(2))==x),peakChannel(2),15];

% plot(x,fitfunction(beta0,x))
hold on 
w = 1./yerr.^2;
w = ones(size(yerr));
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x,y,@fitfunction,beta0,'weights',w);

% errorbar(X,Y,Yerr,'.','markersize',10)
histogram('BinEdges',[0,X(1:end)],'BinCounts',Y,'EdgeColor','none')

plot(x,fitfunction(beta,x),'-.','linewidth',3)
plot(x,beta(1).*x+beta(2)+beta(3).*exp(-((x-beta(4))./(2*beta(5))).^2),'--','linewidth',1)
plot(x,beta(1).*x+beta(2)+beta(6).*exp(-((x-beta(7))./(2*beta(8))).^2),'--','linewidth',1)
plot(x,beta(1).*x+beta(2),'--','linewidth',1)

legend('Data','Total Fit','Ag peak','C peak','Liniar background')

% plot(x,beta(4).*x+beta(5))
us = CovB;
pValue = 1-chi2cdf(MSE*(length(y)-5),(length(y)-5));
P_Value = 1-chi2cdf(MSE*(length(y)-5),(length(y)-5));

disp([name 'Fitted with p-value: ' num2str(pValue) ' with MSE: ' num2str(MSE)])


c2E = @(x) x.*0.76468+12.393;

peakChannelsFitted = [beta(4),beta(7)];
% peakUns(1) = us(1,1);
ci = nlparci(beta,R,'jacobian',J,'alpha',0.35);
peakUns = [(ci(4,2)-ci(4,1))/2,(ci(7,2)-ci(7,1))/2];
peakUns = [us(4),us(7)];

counts1 = abs(2*sqrt(pi)*beta(3)*beta(5));
countUns1 = sqrt(counts1+...
    (2*sqrt(pi)*beta(3)).^2*us(5,5).^2+...
    (2*sqrt(pi)*beta(5)).^2*us(3,3).^2+...
    abs(2.*(2*sqrt(pi)*beta(3)).*(2*sqrt(pi)*beta(5))*us(3,5)));

counts2 = abs(2*sqrt(pi)*beta(6)*beta(8));
countUns2 = sqrt(counts2+...
    (2*sqrt(pi)*beta(6)).^2*us(8,8).^2+...
    (2*sqrt(pi)*beta(8)).^2*us(6,6).^2+...
    abs(2.*(2*sqrt(pi)*beta(6)).*(2*sqrt(pi)*beta(8))*us(6,8)));



data1 = [peakChannelsFitted(1);peakUns(1);[counts1];[countUns1]];
data2 = [peakChannelsFitted(2);peakUns(2);[counts2];[countUns2]];
end


function y = fitfunction(beta,x)
    y = beta(1).*x+beta(2)+beta(3).*exp(-((x-beta(4))./(2.*beta(5))).^2)+beta(6).*exp(-((x-beta(7))./(2.*beta(8))).^2);
    y =y';
end