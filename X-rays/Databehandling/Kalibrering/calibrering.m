clear all; close all;clc;
%% Resulltat 
a = 0.031179;
aUs = 3.605e-06;
b = -0.08267;
bUs = 0.0025276;
coVarAB = -8.2792e-09;
c2E = @(c) a*c+b;
c2EUs = @(c,cUs) sqrt((cUs*a).^2+(c*aUs).^2+bUs.^2+coVarAB.*c);
E2c = @(E) (E-b)/a


%%

filenames = {'Am241.txt','Cs137.txt','Fe55.txt'};
names = {'Am241','Cs137','Fe55'}

% In keV
% peakValues = {[26.350,59.5409],[],[5.8943,6.490]};
% peakValuesUs = [0.001, 0.0001,0.0001,0.0001,0.00001]
% dPeakBorder = [998,1061]
% dPeakValue = 32.1936
% peakBorders = {[835,865;1881,1931],[],[181,202;203,215]};

peakValues = {[13.9,16.84,17.7502,20.7848,26.350,59.5409],[36.6],[5.8943,6.490]};
peakValuesUs = [0.1,0.01,  0.0001, 0.0001, 0.001, 0.0001,   0.1,0.0001,0.0001,0.00001]
dPeakBorder = [998,1061]
dPeakValue = 32.1936
peakBorders = {[430,465;528,559;559,588;661,678;835,865;1881,1931],[1152,1186],[181,202;203,215]};



MSECount = [];
pValue = [];
peakChannel = [];
peakUns = [];
peakE = [];

colorPlot = {'#0072BD','#D95319','#EDB120'}	
	

for i = 1:length(filenames)
    dat = load(['..\..\Kalibrering\' filenames{i}]);
    X = dat(:,1);
    Y = dat(:,2);
    YUs = sqrt(Y)+(Y==0);
        
    % Specifik for spectrum 
    name = names{i};
    peakValue = peakValues{i};
    peakBorder = peakBorders{i};
    
    
    n = length(peakValue);
    
    %%  Passer kalibrerigen?
    figure(i*2-1)
    errorbar(c2E(X),Y,YUs,'.')
    xlabel('Energy (E) [keV]')
    ylabel('Counts (n)')
    set(gca,'FontSize',15) 
    hold on
    for j = 1:n
        plot([peakValue(j),peakValue(j)]+c2EUs([peakValue(j),peakValue(j)],0),[0,max(Y)])
        plot([peakValue(j),peakValue(j)]-c2EUs([peakValue(j),peakValue(j)],0),[0,max(Y)])
    end
    if i ==2
        plot([dPeakValue,dPeakValue]+c2EUs([dPeakValue,dPeakValue],0),[0,max(Y)])
        plot([dPeakValue,dPeakValue]-c2EUs([dPeakValue,dPeakValue],0),[0,max(Y)])
    end
    
    %% Fit data
%     figure(10)
%     errorbar(X,Y,YUs,'.','markersize',10,'color',colorPlot{i})
%     hold on
%     xlabel('Channel number (Ch)')
%     ylabel('Counts (n)')
%     set(gca,'FontSize',15) 
    
    %
    fig = figure(i*2)
    errorbar(X,Y,YUs,'.','markersize',10)
    hold on
    xlabel('Channel number (Ch)')
    ylabel('Counts (n)')
    set(gca,'FontSize',15) 
    title(name)
    % Enkelt peaks
    ax1 = gca;
    hold(ax1, 'on');
    xlim([min(min(peakBorder))-100,max(max(peakBorder))+100])
    ylim([0,max(Y)*1.2])
    
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
        beta0 = [(xmax+xmin)/2,(xmax-xmin)/3,max(y),0,0];
        % plot(x,fitfunction(beta0,x))
        w = 1./yUs.^2;
        [beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x,y,fitfunction,beta0,'weights',w);

        plt = plot(ax1,x,fitfunction(beta,x),'k--','linewidth',2)
        plot(ax1,x,beta(4).*x+beta(5))
        var = CovB/MSE;

        txt = text(ax1,beta(1),max(y)+2,['\leftarrow' num2str(round(peakValue(j),1)) ' MeV']);
        set(txt,'Rotation',90);
        set(txt,'FontSize',13);
%         xan = [0.1,0.2];
%         yan = [0.1,0.5];
%         a = annotation(plt,'textarrow',xan,yan,'String',[num2str(round(peakValue(j),1)) ' MeV']);
% %         a = annotation('textarrow',beta(1),max(y)+2,'String',[num2str(round(peakValue(j),1)) ' MeV']);
%         set(a,'FontSize',13);

        if (i==1)&&(j==1)
            ax2 = axes('Position',[(0.149+0.756*0.4) (0.1481+0.771*0.4) (0.756*(0.6)) (0.771*(0.6))])
            t = title('Fit of one peak', 'Units', 'normalized', 'Position', [0.5, -0.2, 0])
            box on
            hold on
            errorbar(X,Y,YUs,'.','markersize',10)
            plot(x,fitfunction(beta,x),'k--','linewidth',2)
            plot(x,beta(4).*x+beta(5))
            xlim([min(x)-10,max(x)+10])
            

        end
%         figure(10)
%         plot(x,fitfunction(beta,x),'k--','linewidth',2)
%         txt = text(beta(1),max(y)+2,['\leftarrow' num2str(peakValue(j)) ' MeV']);
%         set(txt,'Rotation',90);
%         set(txt,'FontSize',14);
%         
        % plot([beta(1)+us(1,1),beta(1)+us(1,1)],[0,max(y)])
        % plot([beta(1)-us(1,1),beta(1)-us(1,1)],[0,max(y)])

        MSECount(end+1) = MSE;
        pValue(end+1) = 1-chi2cdf(MSE*(length(y)-5),(length(y)-5));
        peakChannel(end+1) = beta(1,1);
        peakUns(end+1) = sqrt(var(1,1));
        peakE(end+1) = peakValue(j);
        
%         figure(i*2)
    end
    %%Double Peak
    if i ==2
        x = X;
        y = Y;
        yerr = YUs;
        xmin = dPeakBorder(1);
        xmax = dPeakBorder(2);
        
        indexes = (x>xmin)&(x<xmax);
        x = x(indexes);
        y = y(indexes);
        yerr = yerr(indexes);
        
        fitfunction= @(beta,x)  beta(3).*exp(-((x-beta(1))./(beta(2))).^2./2)+beta(6).*exp(-((x-beta(4))./(beta(5))).^2./2)+beta(7).*x+beta(8);
        beta0 = [1035,5,146,1023,5,98,0,2];
%         plot(x,fitfunction(beta0,x))
        w = 1./yerr.^2;
        [beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x,y,fitfunction,beta0,'weights',w);

        plot(x,fitfunction(beta,x),'k--','linewidth',2)
        plot(x,beta(7).*x+beta(8)+beta(3).*exp(-((x-beta(1))./(beta(2))).^2./2))
        plot(x,beta(7).*x+beta(8)+beta(6).*exp(-((x-beta(4))./(beta(5))).^2./2))        
        plot(x,beta(7).*x+beta(8))
        var = CovB/MSE;

        txt = text(beta(1),max(y)+2,['\leftarrow' num2str(dPeakValue(j)) ' MeV']);
        set(txt,'Rotation',90);
        set(txt,'FontSize',14);
        
%         figure(10)
%         plot(x,fitfunction(beta,x),'k--','linewidth',2)
%         txt = text(beta(1),max(y)+2,['\leftarrow' num2str(dPeakValue(j)) ' MeV']);
%         set(txt,'Rotation',90);
%         set(txt,'FontSize',14);
%         
        
        MSECount(end+1) = MSE;
        pValue(end+1) = 1-chi2cdf(MSE*(length(y)-8),(length(y)-8));
        peakChannel(end+1) = beta(1);
        peakUns(end+1) = sqrt(var(1,1));
        peakE(end+1) = dPeakValue;
    end
    
    
    
end

x  = peakE;
y = peakChannel;
yUs = peakUns;
yUs = sqrt(peakUns.^2+peakValuesUs.^2./a.^2)

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


