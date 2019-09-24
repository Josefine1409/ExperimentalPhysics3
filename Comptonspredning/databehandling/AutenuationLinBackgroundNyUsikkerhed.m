clear all; close all;clc;
folderName = {'Aluminum','Lead','Cupper','Brass'};
name = {'Al','Pb','Cu','Brass'};
rho = [2.699e+00,1.135E+01,8.960e+00]
mu_rho = [7.5055e-02,11.3663e-2,7.3103E-02]
% mu_rho = [7.4642e-02,11.3663e-2,7.3103E-02]

mu = mu_rho.*rho

channelInterval= [[450,630],[450,630],[450,630],[450,630]];
deltaThiknessLead = [1.13,1.10,1.13,1.01,1.04,1.05,1.05,1.19,1.06,1.21];
plateThiknessLead(1) = 0;
for i=1:length(deltaThiknessLead)
    plateThiknessLead(i+1) = sum(deltaThiknessLead(end-i+1:end));
end
plateThikness = {[0,10.11,20.22,30.34,40.46,50.55,60.65,70.79,80.9]./10,[plateThiknessLead]./10,[0,1.02,2.06,3.05,4.04,5.06,6.10,7.09,8.12,9.11]./10,[0,14.19,29.87,44.42]./10,[0,7.98,16.03,24.02]./10};
plateThiknessUncertanty = 0.02;

for i = 1:length(name)
    counts=zeros(size(plateThikness{i}));
    countUs=zeros(size(plateThikness{i}));
    pathStart = ['..\data\AttenuationCoefficient\' folderName{i} '\AttenuationCoefficient_' name{i} '_'];
        for j = 1:length(plateThikness{i})
            path  = [pathStart num2str(j-1) 'Plates_ch001.txt'];
            [X,Y,Yerr] = hisFraData(path);
            [counts(j),countUs(j)]  = gaussCounter(X,Y,Yerr,400,700,[name{i},': measurement num: ',num2str(j),', T =',num2str(plateThikness{i}(j))]);
%             counts(j) = sum(Y((channelInterval(1,1):channelInterval(1,2))));
%             figure
%             title([name{i},': measurement num: ',num2str(j),', T =',num2str(plateThikness{i}(j))])
%             hold on
%             plot(X,Y,'.')
%             xlabel('Channel')
%             ylabel('Counts')
%             plot(X([channelInterval(1,1):channelInterval(1,2)]),Y([channelInterval(1,1):channelInterval(1,2)]),'*')
        end 
%     figure
%     hej = log(counts./counts(1));
%     plot(plateThikness{i},hej,'.')
    fig = figure
    hold on
    xlabel('Thickness x [cm]')
    ylabel('Counts under peak N')
    set(gca, 'YScale', 'log')
    set(gca,'FontSize',15) 

    x = plateThikness{i};
    y = counts;
    yerr = countUs;
    beta0 = [y(1)-y(end),-0.1];
    linFun =@(beta,x) beta(1).*exp(x.*beta(2));
%     plot(x,linFun(beta0,x))
    w = 1./yerr.^2;
    [beta,R,J,CovB,MSE]  = nlinfit(x,y,@(beta,x) linFun(beta,x),beta0,'weights',w);
%     beta(2)
%     CovB(2,2)
%     MSE
    yerr = sqrt(yerr.^2 + (beta(1).*beta(2).*exp(beta(2).*x)*plateThiknessUncertanty).^2);
    yerr(1) = sqrt(y(1));
    w = 1./yerr.^2;
    [beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x,y,@(beta,x) linFun(beta,x),beta,'weights',w);
%     beta(2)
%         CovB(2,2)
%         MSE

    hold on
    xs = linspace(x(1),x(end),1000);
    plot(xs,linFun(beta,xs),'b','Linewidth',1.5)
    
    %     plot(xs,linFun([beta(1),-mu(i)],xs),'--','Linewidth',1.5)
    
    errorbar(x,y,yerr,'.','markersize',10)
    [Ypred,delta] = nlpredci(@(beta,x) linFun(beta,x),xs,beta,R,'jacobian',J,'alpha',0.35);
    plot(xs,Ypred+delta,'k--','Linewidth',0.5)
    plot(xs,Ypred-delta,'k--','Linewidth',0.5)
    
    legend('Fit','Data','Fit confidence')

    ylim([min(y)+(min(y)-max(y))/0.2-max(yerr),max(y)-(min(y)-max(y))/5])
    us = CovB/MSE;
    US{i} = us;
    mse{i} = MSE;
    BETA{i} = beta;
    pValue= 1-chi2cdf(MSE*(length(y)-2),(length(y)-2));
    PvALUE{i} = pValue;
%     title(['Fit of autenuation for ',name{i}])
    ci = nlparci(beta,R,'jacobian',J,'alpha',0.35)
    
    disp('_____________________________________________________________________________')
    disp(['Fit of atenuation for ',name{i}, ' gives mu = ', num2str(beta(2)),'+-',num2str((ci(2,2)-ci(2,1))/2),' cm^-1'])
%     disp(['and bacground = ', num2str(beta(3)),'+-',num2str(us(3,3)),' counts.'])
    disp(['With MSE = ', num2str(MSE),' with p-value: ',num2str(pValue)])
    saveas(fig,['..\figure\Atenuation',name{i},'.png'])

end



function [X,Y,Yerr] = hisFraData(filename)
delimiter = ' ';
startRow = 6;
formatSpec = '%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
timestamp = dataArray{:, 1};
channel = dataArray{:, 2};
VarName5 = dataArray{:, 3};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;
X = 1:max(channel);
for i = X
    Y(i) = sum(channel==i);
end
Yerr = sqrt(Y) +(Y==0);
end

function [counts,countUns]=gaussCounter(X,Y,Yerr,xmin,xmax,titleName)
figure
hold on
xlabel('Channel Ch ')
ylabel('Counts n')
set(gca,'FontSize',15) 

    
x = X;
y = Y;
yerr = Yerr;
higherIndex = x>xmin;
x = x(higherIndex);
y = y(higherIndex);
yerr = yerr(higherIndex);

lowerIndex = x<xmax;
x = x(lowerIndex);
y = y(lowerIndex);
yerr = yerr(lowerIndex);

beta0 = [(xmin+xmax)/2,(xmax-xmin)/3,max(y),0,20];
% plot(x,fitfunction(beta0,x))
w = 1./yerr.^2;
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x,y,@fitfunction,beta0,'weights',w);
errorbar(X,Y,Yerr,'.','markersize',10)
plot(x,fitfunction(beta,x),'Linewidth',1.5);
plot(x,beta(4).*x+beta(5),'r--','Linewidth',1.5);

[Ypred,delta] = nlpredci(@fitfunction,x,beta,R,'jacobian',J,'alpha',0.35);

plot(x,Ypred+delta,'k--','Linewidth',0.5)
legend('Data','Gauss fit','Background fit','Confidence interval')

plot(x,Ypred-delta,'k--','Linewidth',0.5)


us = CovB/MSE;
us;
pValue = 1-chi2cdf(MSE*(length(y)-5),(length(y)-5));


ci = nlparci(beta,R,'jacobian',J,'alpha',0.35);
us2 = (us(2,2)-us(2,1))/2;
us3 = (us(3,2)-us(3,1))/2;
% background = (xmax-xmin)*beta(5) + 1/2.*beta(4).*(xmax^2-xmin^2)
% title([titleName 'med p-value: ' num2str(pValue)])
counts = abs(2*sqrt(pi)*beta(2)*beta(3));
countUns = sqrt(counts+...
    (2*sqrt(pi)*beta(3)).^2*us2.^2+...
    (2*sqrt(pi)*beta(2)).^2*us3.^2+...
    abs(2.*(2*sqrt(pi)*beta(2)).*(2*sqrt(pi)*beta(3))*us(2,3)));
end
function y = fitfunction(beta,x)
    y = beta(4).*x+beta(5)+abs(beta(3)).*exp(-((x-beta(1))./(2.*abs(beta(2)))).^2);
end