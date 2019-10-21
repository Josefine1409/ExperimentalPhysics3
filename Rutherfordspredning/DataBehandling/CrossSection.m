clear all; close all;clc;
    
Ein = 349.5

c2E = @(x) x*0.76468+12.393
theta =       [30,            40,           50,             60,            70,            75,            110,           120,           130,           140,           150,           160];
peakValues =  {c2E([434,445]),c2E([434,431]),c2E([433,427]),c2E([433,410]),c2E([433,396]),c2E([428,388]),c2E([432,350]),c2E([433,340]),c2E([433,330]),c2E([433,323]),c2E([433,318]),c2E([433,314])};
peakBorders = {[390,480],     [380,490],     [380,480],     [350,460],     [340,460],     [340,455],     [310,455],     [300,455],     [280,455],     [280,455],     [280,455],     [280,455]};

linFun =@(beta,x) (x-beta(2))/beta(1);

dataAg = [];
dataC = [];



for i = 1:length(theta)
    [X,Y,Yerr] = hisFraData(['..\Data\AngularDependency\' num2str(theta(i)) 'degree.asc']);
    [Data1,Data2] = fitGaussInSpectrum(X,Y,Yerr,[num2str(theta(i)) ' degree'],peakValues{i},peakBorders{i});
    dataC = [dataC,Data2];
    dataAg = [dataAg,Data1];

end
%%
EnergyC = c2E(dataC(1,:));
EnergyUnsC = c2E(dataC(2,:));
CountsC = c2E(dataC(3,:));
CountsUnsC = c2E(dataC(4,:));
EnergyAg = c2E(dataAg(1,:));
EnergyUnsAg = c2E(dataAg(2,:));
CountsAg = c2E(dataAg(3,:));
CountsUnsAg = c2E(dataAg(4,:));


%% Maltekode

thetas=[30:10:70 75 110:10:160]./180.*pi;
ts=[300 300 300 300 350 600 600 600 600 600 600 600];
FCs=[46449 51781 60300 65892 81962 34228 103355 80585 101501 106019 102580 37053];

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


GCs=GCs([3:11])
CCs=CCs([3:1:11])
sigmaGCs = sigmaGCs([3:1:11])
sigmaCCs = sigmaCCs([3:1:11])


cs=@(theta,m,E) 1./(sin(theta/2)).^2.*1./(K2(theta,m).*E).^2;
cs=@(theta,m,E) 1./(sin(theta/2)).^4;
% cs=@(theta,m,E) 1./(sin(theta/2)).^2.*1./(E).^2;


% [betaG,RG,JG,CovB,MSE]=nlinfit(thetas,GCs,@(C,thetas) C(1)*cs(thetas,mG,Ein)+C(2),[1,0]);
% % [betaC,RC,JC,CovB,MSE]=nlinfit(thetas,CCs,@(C,thetas) C(1)*cs(thetas,mC,Ein)+C(2),[1,0]);
% [betaG,RG,JG,CovB,MSE]=nlinfit(thetas([1:5,7:11]),GCs,@(C,thetas) C(1)*cs(thetas,mG,Ein)+C(2),[1,0],'weights',1./sigmaGCs.^2);
% [betaC,RC,JC,CovB,MSE]=nlinfit(thetas([1:5,7:11]),CCs,@(C,thetas) C(1)*cs(thetas,mC,Ein)+C(2),[1,0],'weights',1./sigmaCCs.^2);

[betaG,RG,JG,CovB,MSEG]=nlinfit(thetas([3:1:11]),GCs,@(C,thetas) C(2)+ C(1)*cs(thetas,mG,Ein),[1,0]);
[betaC,RC,JC,CovB,MSEC]=nlinfit(thetas([3:1:11]),CCs,@(C,thetas)C(2)+ C(1)*cs(thetas,mC,Ein),[1,0]);

Chi2TestGold = sum((GCs-betaG(1).*cs(thetas([3:1:11]),mG,Ein)-betaG(2)).^2./(sigmaGCs.^2)) % Chi i anden test
nu = length(GCs)-1 %Antal frihedsgrader
Signifikansniveau = 1-chi2cdf(Chi2TestGold,nu) % P værdi
TestN = Chi2TestGold/nu %Ca. test der burde give 1
disp('---------------------')


Chi2TestCarbon = sum((CCs-betaC(1).*cs(thetas([3:1:11]),mC,Ein)-betaC(2)).^2./(sigmaCCs.^2)) % Chi i anden test
nu = length(GCs)-1 %Antal frihedsgrader
Signifikansniveau = 1-chi2cdf(Chi2TestCarbon,nu) % P værdi
TestN = Chi2TestCarbon/nu %Ca. test der burde give 1
disp('---------------------')


Thetas=linspace(thetas(3),thetas(end),1000);

fig = figure
errorbar(theta(3:end),GEs(3:end),sigmaGEs(3:end),'r.','markersize',10)
hold on
xlabel('\theta [deg]')
ylabel('Scattering Energy [KeV]')
title('Scattering on Gold')
plot(Thetas/pi*180,Ein*K2(Thetas,mG),'linewidth',2)
set(gca,'FontSize',15) 
xlim([theta(1),theta(end)])

fig= figure
errorbar(theta(3:end),CEs(3:end),sigmaCEs(3:end),'r.','markersize',10)
hold on
xlabel('\theta [deg]')
ylabel('Scattering Energy [KeV]')
title('Scattering on Carbon')
plot(Thetas/pi*180,Ein*K2(Thetas,mC),'linewidth',2)
set(gca,'FontSize',15) 
xlim([theta(1),theta(end)])

fig = figure
[Ypred,deltaY] = nlpredci(@(C,th) C(2)+C(1)*cs(th,mG,Ein),Thetas,betaG,RG,'jacobian',JG,'alpha',0.35);

hold on
xlabel('\theta [Deg]')
ylabel('Normed Counts')
title('Scattering on Gold')
plot(Thetas/pi*180,Ypred,'linewidth',2)
errorbar(theta([3:1:11]),GCs,sigmaGCs,'r.','markersize',10)

plot(Thetas/pi*180,Ypred+deltaY,'k--','linewidth',1)
plot(Thetas/pi*180,Ypred-deltaY,'k--','linewidth',1)
legend('Fit','Data','Fit confidence')
set(gca,'FontSize',15) 
xlim([theta(1),theta(end)])

fig = figure
[Ypred,deltaY] = nlpredci(@(C,th)C(2)+ C(1)*cs(th,mC,Ein),Thetas,betaC,RC,'jacobian',JC,'alpha',0.35);

hold on
xlabel('\theta [deg]')
ylabel('Normed Counts')
title('Scattering on Carbon')
plot(Thetas/pi*180,Ypred,'linewidth',2)
errorbar(theta([3:1:11]),CCs,sigmaCCs,'r.','markersize',10)

plot(Thetas/pi*180,Ypred+deltaY,'k--','linewidth',1)
plot(Thetas/pi*180,Ypred-deltaY,'k--','linewidth',1)
set(gca,'FontSize',15) 
xlim([theta(1),theta(end)])


legend('Fit','Data','Fit confidence')

%%

fig = figure
ThetasLinG = 1./(sin(Thetas/2)).^4;
thetasLinG = 1./(sin(theta([3:1:11])/180*pi/2)).^4;

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


figure
ThetasLinC = 1./(sin(Thetas/2)).^4;
thetasLinC = 1./(sin(theta([3:1:11])/180*pi/2)).^4;
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

%%




function [X,Y,Yerr] = hisFraData(filename)
% addpath('..\data\Kalibrering')
% delimiter = ' ';
% startRow = 6;
% formatSpec = '%f%f%f%[^\n\r]';
% fileID = fopen(filename,'r');
% dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
% fclose(fileID);
% timestamp = dataArray{:, 1};
% channel = dataArray{:, 2};
% VarName5 = dataArray{:, 3};
% clearvars filename delimiter startRow formatSpec fileID dataArray ans;
% 
% X = 1:max(channel);
% for i = X
%     Y(i) = sum(channel==i);
% end
% Yerr = sqrt(Y) +(Y==0);
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

% 
% figure
% errorbar(X*0.76535+15.78,Y,Yerr,'.')
% % errorbar(X,Y,Yerr,'.','markersize',10)
% 
% xlabel('Energy (E) [MeV]')
% % xlabel('Channel number')
% ylabel('Counts (n)')
% set(gca,'FontSize',15) 
% 
% hold on
% for i = 1:n
%     plot([peakValue(i),peakValue(i)],[0,500])
% end


figure
errorbar(X,Y,Yerr,'.','markersize',10)
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
yerr = Yerr;
higherIndex = x>x2;
x = x(higherIndex);
y = y(higherIndex);
yerr = yerr(higherIndex);

lowerIndex = x<x3;
x = x(lowerIndex);
y = y(lowerIndex);
yerr = yerr(lowerIndex);

peakChannel = (peakValue-15.78)./0.76535;


beta0 = [0,0,y(round(peakChannel(1))==x),peakChannel(1),5,y(round(peakChannel(2))==x),peakChannel(2),15];

for i = 1:n
    plot([peakChannel(i),peakChannel(i)],[0,max(y)])
end

% plot(x,fitfunction(beta0,x))
hold on 
w = 1./yerr.^2;
w = ones(size(yerr));
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x,y,@fitfunction,beta0,'weights',w);
beta(1);
beta(2);
plot(x,fitfunction(beta,x),'linewidth',2)
plot(x,beta(1).*x+beta(2)+beta(3).*exp(-((x-beta(4))./(2*beta(5))).^2),'--')
plot(x,beta(1).*x+beta(2)+beta(6).*exp(-((x-beta(7))./(2*beta(8))).^2),'--')
plot(x,beta(1).*x+beta(2),'--')

% plot(x,beta(4).*x+beta(5))
us = CovB/MSE;
pValue = 1-chi2cdf(MSE*(length(y)-5),(length(y)-5));
P_Value = 1-chi2cdf(MSE*(length(y)-5),(length(y)-5));

disp([name 'Fitted with p-value: ' num2str(pValue) ' with MSE: ' num2str(MSE)])



txt = text(beta(4),y(round(beta(4))==x)+10,['\leftarrow' num2str(beta(4)) '']);
set(txt,'Rotation',90);
set(txt,'FontSize',12);

txt = text(beta(7),y(round(beta(7))==x)+10,['\leftarrow' num2str(beta(7)) '']);
set(txt,'Rotation',90);
set(txt,'FontSize',12);

% plot([beta(1)+us(1,1),beta(1)+us(1,1)],[0,max(y)])
% plot([beta(1)-us(1,1),beta(1)-us(1,1)],[0,max(y)])

peakChannelsFitted = [beta(4),beta(7)];
% peakUns(1) = us(1,1);
ci = nlparci(beta,R,'jacobian',J,'alpha',0.35);
peakUns = [(ci(4,2)-ci(4,1))/2,(ci(7,2)-ci(7,1))/2];
% peakUns(i) = us(1,1);

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