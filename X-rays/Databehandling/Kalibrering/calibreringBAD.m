clear all; close all;clc;
%% Resulltat 
a = 0.03098;
aUs = 3.2083e-06;
b = 0.00086048;
bUs = 0.0018268;

c2E = @(c) a*c+b;
c2EUs = @(c,cUs) sqrt((cUs*a).^2+(c*aUs).^2+bUs.^2);

E2c = @(E) (E-b)/a


%%

filenames = {'Am241.txt','Cs137.txt','Fe55.txt'};
% In keV
peakValues = {[13.9,16.84,17.7502,20.7848,26.350,59.5409],[34.9869],[5.8943,6.49045]};
names = {'Am241','Cs137','Fe55'}
peakBorders = {[438,461;533,559;559,588;661,678;835,865;1881,1931],[1152,1186],[181,202;203,215]};

%[34.2789,34.7197][39.2573]

dPeakValues = {[16.84,17.7502],[30.6251,30.9728],[]};
dPeakBorders = {[533,588],[998,1061],[]};

MSECount = [];
pValue = [];
peakChannel = [];
peakUns = [];
peakE = [];

for i = 1:length(filenames)
    dat = load(['..\..\Kalibrering\' filenames{i}])
    X = dat(:,1)
    Y = dat(:,2)
    YUs = sqrt(Y)+(Y==0)
        
    % Specifik for spectrum 
    name = names{i}
    peakValue = peakValues{i}
    peakBorder = peakBorders{i}
    dPeakValue = dPeakValues{i}
    dPeakBorder = dPeakBorders{i}
    
    
    n = length(peakValue);
    m = length(dPeakValue)/2;
    
    %%  Passer kalibrerigen?
    figure
    errorbar(c2E(X),Y,YUs,'.')
    xlabel('Energy (E) [keV]')
    ylabel('Counts (n)')
    set(gca,'FontSize',15) 
    hold on
    for j = 1:n
        plot([peakValue(j),peakValue(j)]+c2EUs([peakValue(j),peakValue(j)],0),[0,max(Y)])
        plot([peakValue(j),peakValue(j)]-c2EUs([peakValue(j),peakValue(j)],0),[0,max(Y)])
    end
    for j = 1:m
        plot([dPeakValue(j),dPeakValue(j)]+c2EUs([dPeakValue(j),dPeakValue(j)],0),[0,max(Y)])
        plot([dPeakValue(j),dPeakValue(j)]-c2EUs([dPeakValue(j),dPeakValue(j)],0),[0,max(Y)])
        plot([dPeakValue(j+1),dPeakValue(j+1)]+c2EUs([dPeakValue(j+1),dPeakValue(j+1)],0),[0,max(Y)])
        plot([dPeakValue(j+1),dPeakValue(j+1)]-c2EUs([dPeakValue(j+1),dPeakValue(j+1)],0),[0,max(Y)])
    end
    
    %% Fit data
    figure
    errorbar(X,Y,YUs,'.','markersize',10)
    hold on
    xlabel('Channel number (Ch)')
    ylabel('Counts (n)')
    set(gca,'FontSize',15) 
    title(name)
    % Enkelt peaks
    for j = 1:n
        %Udvælger data
        x = X;
        y = Y;
        yUs = YUs;
        xmin = peakBorder(j,1);
        xmax = peakBorder(j,2);
        
        indexes = (x>xmin)&(x<xmax);
        x = x(indexes);
        y = y(indexes);
        yUs = yUs(indexes);
        
        fitfunction= @(beta,x)  beta(4).*x+beta(5)+beta(3).*exp(-((x-beta(1))./(beta(2))).^2./2);
        beta0 = [(xmax+xmin)/2,(xmax-xmin)/3,max(y),0,mean(Y(1:(floor(end/3))))];
        % plot(x,fitfunction(beta0,x))
        w = 1./yUs.^2;
        [beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x,y,fitfunction,beta0,'weights',w);

        plot(x,fitfunction(beta,x),'k--','linewidth',2)
        plot(x,beta(4).*x+beta(5))
        var = CovB/MSE;

        txt = text(beta(1),max(y)+2,['\leftarrow' num2str(peakValue(j)) ' MeV']);
        set(txt,'Rotation',90);
        set(txt,'FontSize',14);

        % plot([beta(1)+us(1,1),beta(1)+us(1,1)],[0,max(y)])
        % plot([beta(1)-us(1,1),beta(1)-us(1,1)],[0,max(y)])

        MSECount(end+1) = MSE;
        pValue(end+1) = 1-chi2cdf(MSE*(length(y)-5),(length(y)-5));
        peakChannel(end+1) = beta(1,1);
        peakUns(end+1) = sqrt(var(1,1));
        peakE(end+1) = peakValue(j)
    end
    %%Double Peaks
    for j = 1:m
        %Udvælger data
        x = X;
        y = Y;
        yerr = YUs;
        xmin = dPeakBorder(j,1);
        xmax = dPeakBorder(j,2);
        
        indexes = (x>xmin)&(x<xmax);
        x = x(indexes);
        y = y(indexes);
        yerr = yerr(indexes);
        
        fitfunction= @(beta,x)  beta(3).*exp(-((x-beta(1))./(beta(2))).^2./2)+beta(6).*exp(-((x-beta(4))./(beta(5))).^2./2)+beta(7).*x+beta(8);
        beta0 = [E2c(dPeakValue(j)),(xmax-xmin)/9,y((round(end*1/4))),E2c(dPeakValue(j+1)),(xmax-xmin)/9,y((round(end*3/4))),0,mean(Y(1:(floor(end/3))))];
        plot(x,fitfunction(beta0,x))
        w = 1./yerr.^2;
        [beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x,y,fitfunction,beta0,'weights',w);

%         plot(x,fitfunction(beta,x),'k--','linewidth',2)
%         plot(x,beta(4).*x+beta(5)+beta(3).*exp(-((x-beta(1))./(beta(2))).^2./2))
%         plot(x,beta(4).*x+beta(5)+beta(6).*exp(-((x-beta(4))./(beta(5))).^2./2))        
%         plot(x,beta(4).*x+beta(5))
%         var = CovB/MSE;
% 
%         txt = text(beta(1),max(y)+2,['\leftarrow' num2str(dPeakValue(j)) ' MeV']);
%         set(txt,'Rotation',90);
%         set(txt,'FontSize',14);
%         
%         txt = text(beta(4),max(y)+2,['\leftarrow' num2str(dPeakValue(j+1)) ' MeV']);
%         set(txt,'Rotation',90);
%         set(txt,'FontSize',14);
%         
%         % plot([beta(1)+us(1,1),beta(1)+us(1,1)],[0,max(y)])
%         % plot([beta(1)-us(1,1),beta(1)-us(1,1)],[0,max(y)])
% 
%         MSECount(end+1) = MSE;
%         pValue(end+1) = 1-chi2cdf(MSE*(length(y)-5),(length(y)-5));
%         peakChannel(end+1) = beta(1);
%         peakUns(end+1) = sqrt(var(1,1));
%         peakE(end+1) = dPeakValue(j)
% 
%         peakChannel(end+1) = beta(4);
%         peakUns(end+1) = sqrt(var(4));
%         peakE(end+1) = dPeakValue(j+1)
    end
    
    
    
end

x  = peakE;
y = peakChannel;
yUs = peakUns;
Es = linspace(min(x),max(x),1000);      

beta0 = [1,2];
linFun =@(beta,x) (x-beta(2))/beta(1);
w = 1./yUs.^2;
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x,y,@(beta,x) linFun(beta,x),beta0,'weights',w);

var = CovB/MSE;
[Ypred,delta] = nlpredci(@(beta,x) linFun(beta,x),Es,beta,R,'jacobian',J,'alpha',0.32);

figure
hold on
xlabel('Energy (E) [keV]')
ylabel('Channel number (Ch)')
set(gca,'FontSize',15) 

plot(Es,Ypred,'-','linewidth',1)
errorbar(x,y',yUs,'.','markersize',8)

plot(Es,Ypred+delta,'b--','linewidth',1)
plot(Es,Ypred-delta,'b--','linewidth',1)
legend('Fit','Data','68% Confidence','Location','northwest')

figure
xlabel('Energy (E) [keV]')
ylabel('Residual Channel number (\DeltaCh)')
set(gca,'FontSize',15) 
hold on

plot(Es,0.*Es,'b','linewidth',1)

[YpredMeasurement,deltaM] = nlpredci(@(beta,x) linFun(beta,x),x,beta,R,'jacobian',J,'alpha',0.35);
errorbar(x,y-YpredMeasurement',yUs,'.r','markersize',8)
plot(Es,delta,'k--','linewidth',1)
plot(Es,-delta,'k--','linewidth',1)

legend('Fit','Data','68% confidence','Location','southwest')

var = CovB/MSE;
MSE
pValue = 1-chi2cdf(MSE*(length(y)-2),(length(y)-2))

% disp(['Energi som funktion af channel: E= ' num2str(beta(1)) '+-' num2str((ci(1,2)-ci(1,1))/2) 'Channel+' num2str(beta(2)) '+-' num2str((ci(2,2)-ci(2,1))/2) 'MeV'])
disp(['Energi som funktion af channel: E= ' num2str(beta(1)) '+-' num2str(sqrt(var(1,1))) 'Channel+' num2str(beta(2)) '+-' num2str(sqrt(var(2,2))) 'keV'])
