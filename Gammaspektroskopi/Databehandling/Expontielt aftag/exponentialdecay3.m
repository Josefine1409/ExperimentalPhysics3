clear all; close all;clc;
%Energy calibration 
c2E = @(c) 0.73642*c +0.06;
c2EUs = @(c,cUs) sqrt((cUs*0.73642).^2+(c*0.00004).^2+0.05.^2);

%Efficientcy calibration 
rho = @(E) 1.9.*E.^(-0.95)+0.0010;
rhoUs = @(E) sqrt(...
    0.0002^2+...
    (0.5.*E.^(-0.95)).^2+...
    (-0.95.*1.9.*E.^(-1.95).*0.04).^2);
counts2Decay = @(counts,E) counts./rho(E) ;
countsUs2Decay =@(counts,countsUs,E) sqrt((countsUs./rho(E)).^2+(counts./rho(E).^2.*rhoUs(E)));

%precent of decay
BR = 0.9882;
BRUs = 0;

filenames = {'Mn56_Day1_ch000.txt','Mn56_Day2_ch000.txt','Mn56_Day3_2detectors_ch000.txt','Mn56_Day3_2detectors_ch002.txt'};
startGuess = [1150,1150,1150,935];

%starttime in seconds
startTime = [22.07,25.36,21.2,21.2];
startTimeUs = [1,1,1,1];
%Mass of sample in g
%%
mass = [2.50,2.50,4.10];
massUs = [0.05,0.05,0.05];
%%
%%Constants
M = 54.93804;
NA = 6.02214076*10^23; 
CrossSection = 13.36 *10^-28;
CrossSectionUs = 0.05*10^-28;
DATA = {};
for i = 1:3
   DATA{i} =  load(['..\..\Manganese\' filenames{i}]);
end

%%
close all

for i=1:3
data = load(['..\..\Manganese\' filenames{i}]);

timestamps = data(:,1);
chanelEvent = data(:,2);
n = length(chanelEvent);

% DeltaN = floor(n/30)

DeltaN = 100000;
lowCh  = startGuess(i)-50;
highCh  = startGuess(i)+50;

t = zeros(size(floor(n/DeltaN)));
deltaT = zeros(size(floor(n/DeltaN)));
totDecay = zeros(size(floor(n/DeltaN)));
totDecayUs = zeros(size(floor(n/DeltaN)));
for j = 1:(floor(n/DeltaN))
    indexes = ((j-1).*DeltaN+1):(j.*DeltaN);
    t(j) = (timestamps(indexes(end))+timestamps(indexes(1)))/2;
    deltaT(j) = (timestamps(indexes(end))-timestamps(indexes(1)))/2;
    
    
    c = chanelEvent(indexes);
    [N,edges] = histcounts(c,lowCh:(highCh+1));
    
%     figure(20+j);
    background = (mean([N(1:30)])+mean([N(70:end)]))/2;
    backgroundUs = (std([N(1:30)])+std([N(70:end)]))/2;
%     histogram([N(1:50)])
    
    chN = lowCh:highCh;
    E = c2E(chN);
    Eedges = c2E(edges-1/2);
    decays = counts2Decay(N,E);
    NUs = sqrt(abs(N-background))+backgroundUs;
    decaysUs = countsUs2Decay(N,NUs,E);
    
    %%Fit
    weights = 1./decaysUs.^2;
    [M,I] = max(decays);
    beta0 =[chN(I),5,M,0,background];
    
    fitfunction = @(beta,x) beta(3).*exp(-((x-beta(1))./(beta(2))).^2./2)+beta(4).*x+beta(5);
    [beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(chN,decays,fitfunction,beta0,'weights',weights);
%     ci = nlparci(beta,R,'jacobian',J,'alpha',0.35);
%     sigmaVar = (ci(2,2)-ci(2,1))/2;
%     heightVar = (ci(3,2)-ci(3,1))/2;
%     
    Epeak = c2E(beta(1));
    
    counts = sqrt(2*pi)*beta(3)*beta(2);
%     countsUs = sqrt(...
%         (sqrt(2*pi)*beta(3)*sigmaUs).^2 +...
%         (sqrt(2*pi)*heightUs*beta(2)).^2);

    
    us = CovB/MSE;
    sigmaVar = us(2,2);
    heightVar = us(3,3);
    covarians = us(2,3);
    countsUs = sqrt(...
        (sqrt(2*pi)*beta(3)).^2*sigmaVar +...
        (sqrt(2*pi)*beta(2)).^2*heightVar+...
        4*pi*beta(3)*beta(2)*covarians);
    
%     totDecay(j) = counts2Decay(counts,Epeak);
%     totDecayUs(j) = countsUs2Decay(counts,countsUs,Epeak);
      totDecay(j) = counts;
      totDecayUs(j) = countsUs;
      
      if j==1
        %Plot
        figure(i)
        hold on
        xlabel('E [keV]')
        xlabel('Channel')
        ylabel('Decays')
        set(gca,'FontSize',15)
        errorbar(chN,decays,decaysUs,'.')
%         histogram('BinEdges',Eedges,'BinCounts',decays);
        plot(chN,fitfunction(beta,chN),'k--','linewidth',2);
        plot(chN,beta(4).*chN+beta(5),'k:','linewidth',1);
      end
      if j>1 && (totDecay(j)-totDecay(j-1))>1e5
        %Plot
        figure(j+i*30+100)
        hold on
        xlabel('E [keV]')
        ylabel('Decays')
        set(gca,'FontSize',15)
        errorbar(E,decays,decaysUs,'.')
%         histogram('BinEdges',Eedges,'BinCounts',decays);
        plot(E,fitfunction(beta,chN),'k--','linewidth',2);
        plot(E,beta(4).*chN+beta(5),'k:','linewidth',1);
      end
      

end
%%
%%Omset tid til timer
th = t/(10^8*60*60)+startTime(i)/(60*60);
deltaTh = deltaT/(10^8*60*60);


dDecay_dt = totDecay./deltaTh;
% Hvordan får jeg usikkerhed fra at jeg intigrere over tid?
dDecayUs_dt =totDecayUs./deltaTh;

figure(i+4);
hold on
xlabel('t [h]')
ylabel('Decays pr hour [1/h]')
set(gca,'FontSize',15) 
errorbar(th,dDecay_dt,dDecayUs_dt,'.');

figure(10);
hold on

errorbar(th,dDecay_dt/mass(i),dDecayUs_dt/mass(i),'.');


%fjerner et punkt som svare til påfyldning af detecroen
if(i == 4)
    th = th([1:17 19:end])
    deltaTh = deltaTh([1:17 19:end])
    dDecay_dt = dDecay_dt([1:17 19:end])
    dDecayUs_dt = dDecayUs_dt([1:17 19:end])
end


weights = 1./dDecayUs_dt.^2;
beta0 = [2.5,dDecay_dt(1)];
fitfunction = @(beta,x) beta(2).*exp(-x./beta(1));
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(th,dDecay_dt,fitfunction,beta0,'weights',weights);

BETA(i) = beta(2);
BETA1(i) = beta(1);

ci = nlparci(beta,R,'jacobian',J,'alpha',0.35);
disp(['Maaling: ' num2str(i) '_____________________________________'])
MSE = MSE
P_Value = 1-chi2cdf(MSE*(length(th)-2),(length(th)-2))


us = CovB/MSE;
activety = beta(2);
activetyFitUs = sqrt(us(2,2));

activetyUs(i) = sqrt(activetyFitUs^2+(-activety/beta(1)*startTimeUs(i)/(60*60))^2);
Jh = (activety./BR).*M./(mass(i).*NA.*CrossSection)
JhUs = sqrt(...
    (activetyUs(i)/BR*M/(mass(i)*NA*CrossSection))^2+...
    (activety/BR^2*M/(mass(i)*NA*CrossSection)*BRUs)^2+...
    (activety/BR*M/(mass(i)^2*NA*CrossSection)*massUs(i))^2+...
    (activety/BR*M/(mass(i)*NA*CrossSection^2)*CrossSectionUs)^2);

Jcms = Jh/(60^2*100^2)
JcmsUs = JhUs/(60^2*100^2)

tau = beta(1)
tauUs = (ci(1,2)-ci(1,1))/2

halftime = tau*log(2)
halftimeUs = tauUs*log(2)


ths = linspace(0,max(th),100);
figure(i+4);
plot(ths,fitfunction(beta,ths),'-','linewidth',1)
figure(10);
plot(ths,fitfunction(beta,ths)/mass(i),'-','linewidth',1)


JhList(i) = Jh;
JhUsList(i) = JhUs;


halftimeList(i) = halftime;
halftimeUsList(i) = halftimeUs;

P_ValueList(i) = P_Value;
MSEList(i) =MSE;
end

%%
figure(10)
legend('Data day 1','Fit day 1','Data day 2','Fit day 2','Data day 3','Fit day 3')
set(gca,'yscale','log')
xlabel('t [h]')
ylabel('Decays pr hour pr mass [h^{-1} g^{-1}]')
set(gca,'FontSize',15) 


P_ValueList
MSEList


halftimeList
halftimeUsList

WhalfTime = 1./halftimeUsList(1:3).^2
halftime = halftimeList(1:3)*WhalfTime'./sum(WhalfTime)
halftimeUs = (sqrt(sum(WhalfTime)))^-1
std(halftimeUsList(1:3))


%%
disp('_________________________________')
Jh = (BETA./BR).*M./(mass(1:3).*NA.*CrossSection)
JhUs = sqrt(...
    (activetyUs./BR.*M./(mass(1:3)*NA*CrossSection)).^2+...
    (BETA./BR^2.*M./(mass(1:3)*NA*CrossSection)*BRUs).^2+...
    (BETA./BR.*M./(mass(1:3).^2.*NA*CrossSection).*massUs(1:3)).^2+...
    (BETA./BR.*M./(mass(1:3)*NA*CrossSection^2)*CrossSectionUs).^2);


Jscm = Jh/(60*60*100^2)
JscmUs = JhUs/(60*60*100^2)


WJ = 1./JscmUs(1:3).^2
J = Jscm(1:3)*WJ'./sum(WJ)
JUs = (sqrt(sum(WJ)))^-1
mean(Jscm)
std(Jscm)