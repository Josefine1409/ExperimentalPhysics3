clear all; close all;clc;
    
Ein = 349.5
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
    peakValues{i} = [Ein.*K2(thetas(i),mG),Ein.*K2(thetas(i),mC)];
end

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

ts=[300 300 300 300 350 600 600 600 600 600 600 600];
FCs=[46449 51781 60300 65892 81962 34228 103355 80585 101501 106019 102580 37053];

GCs=CountsAg;
CCs=CountsC;
GEs=EnergyAg;
CEs=EnergyC;

sigmaGCs=CountsUnsAg;
sigmaCCs=CountsUnsC;
sigmaGEs=EnergyUnsAg;
sigmaCEs=EnergyUnsC;

GCs=GCs./FCs;
CCs=CCs./FCs;
sigmaGCs = sigmaGCs./FCs;
sigmaCCs = sigmaCCs./FCs;




cs=@(theta,m,E) 1./(sin(theta/2)).^2.*1./(K2(theta,m).*E).^2;



[betaG,RG,JG,CovB,MSE]=nlinfit(thetas,GCs,@(C,thetas) C(1)*cs(thetas,mG,Ein)+C(2),[1,0]);
[betaC,RC,JC,CovB,MSE]=nlinfit(thetas,CCs,@(C,thetas) C(1)*cs(thetas,mC,Ein)+C(2),[1,0]);



Thetas=linspace(thetas(1),thetas(end),1000);

figure
errorbar(thetas,GEs,sigmaGEs,'r.','markersize',10)
hold on
xlabel('Scattering Angle [deg]')
ylabel('Scattering Energy [KeV]')
title('Scattering on Carbon')
plot(Thetas,Ein*K2(Thetas,mG),'linewidth',2)

figure
errorbar(thetas,CEs,sigmaCEs,'r.','markersize',10)
hold on
xlabel('Scattering Angle [deg]')
ylabel('Scattering Energy [KeV]')
title('Scattering on Carbon')
plot(Thetas,Ein*K2(Thetas,mC),'linewidth',2)


figure
[Ypred,deltaY] = nlpredci(@(C,th) C(1)*cs(th,mG,Ein)+C(2),Thetas,betaG,RG,'jacobian',JG,'alpha',0.35);

errorbar(thetas,GCs,sigmaGCs,'r.','markersize',10)
hold on
xlabel('Scattering Angle')
ylabel('Normed Counts')
title('Scattering on Gold')
plot(Thetas,Ypred,'linewidth',2)

plot(Thetas,Ypred+deltaY,'k--','linewidth',1)
plot(Thetas,Ypred-deltaY,'k--','linewidth',1)
legend('Fit','Data','Fit confidence','Location','southwest')


figure
[Ypred,deltaY] = nlpredci(@(C,th) C(1)*cs(th,mC,Ein)+C(2),Thetas,betaC,RC,'jacobian',JC,'alpha',0.35);

errorbar(thetas,CCs,sigmaCCs,'r.','markersize',10)
hold on
xlabel('Scattering Angle [deg]')
ylabel('Normed Counts')
title('Scattering on Carbon')
plot(Thetas,Ypred,'linewidth',2)

plot(Thetas,Ypred+deltaY,'k--','linewidth',1)
plot(Thetas,Ypred-deltaY,'k--','linewidth',1)
legend('Fit','Data','Fit confidence','Location','southwest')




%%




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

yerr = Yerr+(Yerr==0);


higherIndex = x>x2;
x = x(higherIndex);
y = y(higherIndex);
yerr = yerr(higherIndex);

lowerIndex = x<x3;
x = x(lowerIndex);
y = y(lowerIndex);
yerr = yerr(lowerIndex);
% c2E = @(x) x.*0.76468+12.393

peakChannel = (peakValue-12.393)./0.76468;


beta0 = [0,0,y(round(peakChannel(1))==x),5,y(round(peakChannel(2))==x),15,peakChannel(2)];

for i = 1:n
    plot([peakChannel(i),peakChannel(i)],[0,max(y)])
end

% plot(x,fitfunction(beta0,x))
hold on 
w = 1./yerr.^2;
w = ones(size(yerr));

fitFun = @(beta,x) fitfunction(beta,peakChannel(1),x) 
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x,y,fitFun,beta0,'weights',w);
beta(1);
beta(2);
plot(x,fitFun(beta,x),'linewidth',2)
plot(x,beta(1).*x+beta(2)+beta(3).*exp(-((x-peakChannel(1))./(2*beta(4))).^2),'--')
plot(x,beta(1).*x+beta(2)+beta(5).*exp(-((x-peakChannel(2))./(2*beta(6))).^2),'--')
plot(x,beta(1).*x+beta(2),'--')

% plot(x,beta(4).*x+beta(5))
us = CovB/MSE;
pValue = 1-chi2cdf(MSE*(length(y)-5),(length(y)-5));
P_Value = 1-chi2cdf(MSE*(length(y)-5),(length(y)-5));

disp([name 'Fitted with p-value: ' num2str(pValue) ' with MSE: ' num2str(MSE)])



txt = text(peakChannel(1),y(round(peakChannel(1))==x)+10,['\leftarrow' num2str(peakChannel(1)) '']);
set(txt,'Rotation',90);
set(txt,'FontSize',12);

txt = text(peakChannel(1),y(round(peakChannel(1))==x)+10,['\leftarrow' num2str(peakChannel(1)) '']);
set(txt,'Rotation',90);
set(txt,'FontSize',12);

% plot([beta(1)+us(1,1),beta(1)+us(1,1)],[0,max(y)])
% plot([beta(1)-us(1,1),beta(1)-us(1,1)],[0,max(y)])

peakChannelsFitted = [peakChannel(1),beta(7)];
% peakUns(1) = us(1,1);
ci = nlparci(beta,R,'jacobian',J,'alpha',0.35);
% peakUns = [(ci(4,2)-ci(4,1))/2,(ci(7,2)-ci(7,1))/2];
% peakUns(i) = us(1,1);

counts1 = abs(2*sqrt(pi)*beta(3)*beta(4));
countUns1 = sqrt(counts1+...
    (2*sqrt(pi)*beta(3)).^2*us(4,4).^2+...
    (2*sqrt(pi)*beta(4)).^2*us(3,3).^2+...
    abs(2.*(2*sqrt(pi)*beta(3)).*(2*sqrt(pi)*beta(4))*us(3,4)));

counts2 = abs(2*sqrt(pi)*beta(5)*beta(6));
countUns2 = sqrt(counts2+...
    (2*sqrt(pi)*beta(5)).^2*us(6,6).^2+...
    (2*sqrt(pi)*beta(6)).^2*us(5,5).^2+...
    abs(2.*(2*sqrt(pi)*beta(5)).*(2*sqrt(pi)*beta(6))*us(5,6)));



data1 = [peakChannelsFitted(1);1;[counts1];[countUns1]];
data2 = [peakChannelsFitted(2);CovB(7,7);[counts2];[countUns2]];
end


function y = fitfunction(beta,p1,x)
    y = beta(3).*exp(-((x-p1)./(2.*beta(4))).^2)+beta(5).*exp(-((x-beta(7))./(2.*beta(6))).^2);
    y =y';
end