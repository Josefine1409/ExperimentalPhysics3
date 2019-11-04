clear all; close all;clc;
    
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

for i = 1:1
    figure
    [X,Y,Yerr] = hisFraData(['..\Data\AngularDependency\' num2str(theta(i)) 'degree.asc']);
    W = 1/Yerr.^2;
    peakValue = peakValues{i}
    peakChannel = (peakValues{i}-11.5073)./0.76675;
    pBorders = peakBorders{i}
    Xs = pBorders(1):pBorders(2)
    Ys = Y(Xs)
    
    a1 =    6.58e+04  (4.044e+04, 9.117e+04)
    a2 =   3.086e+04  (5640, 5.609e+04)
    c1 =       3.459  (-5.533, 12.45)
    c2 =    0.001545  (-3904, 3904)
    sigma1 =       8.793  (8.074, 9.512)
    sigma2 =       5.752  (4.536, 6.968)
    x01 =       435.4  (435, 435.8)
    x02 =       434.2  (433.7, 434.8)
    
    pb = peakBorders{i}
    start = [peakChannel(1),14.4,peakChannel(2),16.8]
%     (signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,autozero,fixedparameters,plots,bipolar,minwidth,DELTA,clipheight)
%     [FitResults,GOF,baseline,coeff,residual,xi,yi,bootstrap]= peakfit([X;Y'],(pb(1)+pb(2))/2,(-pb(1)+pb(2)),2,1,'',2,start,1,[],1,0,7);
%     FitResults
%     
%     EnergyC(i) = c2E(FitResults(2,2));
%     EnergyUnsC(i) = (bootstrap(2,2+4))*0.76675;
%     CountsC(i) = (FitResults(2,4));
%     CountsUnsC(i) = (bootstrap(2,4+4));
%     EnergyAg(i) = c2E(FitResults(1,2));
%     EnergyUnsAg(i) = (bootstrap(1,2+4))*0.76675;
%     CountsAg(i) = (FitResults(1,4));
%     CountsUnsAg(i) = (bootstrap(1,4+4));
  
end

%% Maltekode

% ts=[300 300 300 300 350 600 600 600 600 600 600 600];
% FCs=[46449 51781 60300 65892 81962 34228 103355 80585 101501 106019 102580 37053];
% 
% GCs=CountsAg;
% CCs=CountsC;
% GEs=EnergyAg;
% CEs=EnergyC;
% 
% sigmaGCs=CountsUnsAg;
% sigmaCCs=CountsUnsC;
% sigmaGEs=EnergyUnsAg;
% sigmaCEs=EnergyUnsC;
% 
% GCs=GCs./FCs;
% CCs=CCs./FCs;
% sigmaGCs = sigmaGCs./FCs;
% sigmaCCs = sigmaCCs./FCs;
% 
% GCs=GCs([1:5,7:11])
% CCs=CCs([1:5,7:11])
% sigmaGCs = sigmaGCs([1:5,7:11])
% sigmaCCs = sigmaCCs([1:5,7:11])
% 
% 
% cs=@(theta,m,E) 1./(sin(theta/2)).^4.*1./(K2(theta,m).*E).^2;
% 
% % [betaG,RG,JG,CovB,MSE]=nlinfit(thetas,GCs,@(C,thetas) C(1)*cs(thetas,mG,Ein)+C(2),[1,0]);
% % [betaC,RC,JC,CovB,MSE]=nlinfit(thetas,CCs,@(C,thetas) C(1)*cs(thetas,mC,Ein)+C(2),[1,0]);
% [betaG,RG,JG,CovB,MSE]=nlinfit(thetas([1:5,7:11]),GCs,@(C,thetas) C(1)*cs(thetas,mG,Ein)+C(2),[1,0]);
% [betaC,RC,JC,CovB,MSE]=nlinfit(thetas([1:5,7:11]),CCs,@(C,thetas) C(1)*cs(thetas,mC,Ein)+C(2),[1,0]);
% 
% Thetas=linspace(thetas(1),thetas(end),1000);
% 
% figure
% errorbar(theta,GEs,sigmaGEs,'r.','markersize',10)
% hold on
% xlabel('Scattering Angle [deg]')
% ylabel('Scattering Energy [KeV]')
% title('Scattering on Gold')
% plot(Thetas/pi*180,EoutNy(Ein,Thetas,1),'linewidth',2)
% plot(Thetas/pi*180,Ein*K2(Thetas,mG),'linewidth',2)
% 
% figure
% errorbar(theta,CEs,sigmaCEs,'r.','markersize',10)
% hold on
% xlabel('Scattering Angle [deg]')
% ylabel('Scattering Energy [KeV]')
% title('Scattering on Carbon')
% plot(Thetas/pi*180,EoutNy(Ein,Thetas,2),'linewidth',2)
% plot(Thetas/pi*180,Ein*K2(Thetas,mC),'linewidth',2)
% 
% % 
% % figure
% % [Ypred,deltaY] = nlpredci(@(C,th) C(1)*cs(th,mG,Ein)+C(2),Thetas,betaG,RG,'jacobian',JG,'alpha',0.35);
% % 
% % errorbar(thetas,GCs,sigmaGCs,'r.','markersize',10)
% % hold on
% % xlabel('Scattering Angle')
% % ylabel('Normed Counts')
% % title('Scattering on Gold')
% % plot(Thetas,Ypred,'linewidth',2)
% % 
% % plot(Thetas,Ypred+deltaY,'k--','linewidth',1)
% % plot(Thetas,Ypred-deltaY,'k--','linewidth',1)
% % legend('Fit','Data','Fit confidence','Location','southwest')
% % 
% % 
% % figure
% % [Ypred,deltaY] = nlpredci(@(C,th) C(1)*cs(th,mC,Ein)+C(2),Thetas,betaC,RC,'jacobian',JC,'alpha',0.35);
% % 
% % errorbar(thetas,CCs,sigmaCCs,'r.','markersize',10)
% % hold on
% % xlabel('Scattering Angle [deg]')
% % ylabel('Normed Counts')
% % title('Scattering on Carbon')
% % plot(Thetas,Ypred,'linewidth',2)
% % 
% % plot(Thetas,Ypred+deltaY,'k--','linewidth',1)
% % plot(Thetas,Ypred-deltaY,'k--','linewidth',1)
% % legend('Fit','Data','Fit confidence','Location','southwest')
% 
% figure
% [Ypred,deltaY] = nlpredci(@(C,th) C(1)*cs(th,mG,Ein)+C(2),Thetas,betaG,RG,'jacobian',JG,'alpha',0.35);
% 
% hold on
% xlabel('Scattering Angle')
% ylabel('Normed Counts')
% title('Scattering on Gold')
% plot(Thetas/pi*180,Ypred,'linewidth',2)
% errorbar(theta([1:5,7:11]),GCs,sigmaGCs,'r.','markersize',10)
% 
% plot(Thetas/pi*180,Ypred+deltaY,'k--','linewidth',1)
% plot(Thetas/pi*180,Ypred-deltaY,'k--','linewidth',1)
% legend('Fit','Data','Fit confidence','Location','southwest')
% 
% 
% figure
% [Ypred,deltaY] = nlpredci(@(C,th) C(1)*cs(th,mC,Ein)+C(2),Thetas,betaC,RC,'jacobian',JC,'alpha',0.35);
% 
% hold on
% xlabel('Scattering Angle [deg]')
% ylabel('Normed Counts')
% title('Scattering on Carbon')
% plot(Thetas/pi*180,Ypred,'linewidth',2)
% errorbar(theta([1:5,7:11]),CCs,sigmaCCs,'r.','markersize',10)
% 
% plot(Thetas/pi*180,Ypred+deltaY,'k--','linewidth',1)
% plot(Thetas/pi*180,Ypred-deltaY,'k--','linewidth',1)
% 
% 
% legend('Fit','Data','Fit confidence','Location','southwest')
% 

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
