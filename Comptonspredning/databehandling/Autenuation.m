clear all; close all;clc;
folderName = {'Aluminum','Lead','Cupper','Brass'};
name = {'Al','Pb','Cu','Brass'};
rho = [2.699e+00,1.135E+01,8.960e+00]
mu_rho = [7.5055e-02,11.3663e-2,7.3103E-02]
mu = mu_rho.*rho

channelInterval= [[450,630],[450,630],[450,630],[450,630]];
deltaThiknessLead = [1.13,1.10,1.13,1.01,1.04,1.05,1.05,1.19,1.06,1.21];
plateThiknessLead(1) = 0;
for i=1:length(deltaThiknessLead)
    plateThiknessLead(i+1) = sum(deltaThiknessLead(end-i+1:end));
end
plateThikness = {[0,10.11,20.22,30.34,40.46,50.55,60.65,70.79,80.9]./10,[plateThiknessLead]./10,[0,1.02,2.06,3.05,4.04,5.06,6.10,7.09,8.12,9.11]./10,[0,14.19,29.87,44.42]./10,[0,7.98,16.03,24.02]./10};
plateThiknessUncertanty = 0.01;

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
    figure
    hej = log(counts./counts(1));
    plot(plateThikness{i},hej,'.')
    figure
    hold on
    xlabel('Thikness of material [cm]')
    ylabel('Counts under peak')
    x = plateThikness{i};
    y = counts;
    yerr = countUs;
    beta0 = [y(1)-y(end),-0.1];
    linFun =@(beta,x) beta(1).*exp(x.*beta(2));
%     plot(x,linFun(beta0,x))
    w = 1./yerr.^2;
    [beta,~,~,CovB,MSE]  = nlinfit(x,y,@(beta,x) linFun(beta,x),beta0,'weights',w);
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

    errorbar(x,y,yerr,'.')
    hold on
    xs = linspace(x(1),x(end)+2,1000);
    plot(xs,linFun(beta,xs))
    us = CovB/MSE;
    US{i} = us;
    mse{i} = MSE;
    BETA{i} = beta;
    pValue= 1-chi2cdf(MSE*(length(y)-2),(length(y)-2));
    PvALUE{i} = pValue;
    title(['Fit of autenuation for material ',name{i},', P-value =',num2str(pValue)])
    
    disp('_____________________________________________________________________________')
    disp(['Fit of atenuation for ',name{i}, ' gives mu = ', num2str(beta(2)),'+-',num2str(us(2,2)),' cm^-1'])
%     disp(['and bacground = ', num2str(beta(3)),'+-',num2str(us(3,3)),' counts.'])
    disp(['With MSE = ', num2str(MSE),' with p-value: ',num2str(pValue)])

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
errorbar(X,Y,Yerr,'.')
hold on
xlabel('Channel')
ylabel('Counts')

    
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

beta0 = [(xmin+xmax)/2,(xmax-xmin)/3,max(y),0,20,0];
plot(x,fitfunction(beta0,x))
w = 1./yerr.^2;
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x,y,@fitfunction,beta0,'weights',w);
plot(x,fitfunction(beta,x));
plot(x,beta(4).*x+beta(5)+beta(6).*x.^2);
us = CovB/MSE;
pValue = 1-chi2cdf(MSE*(length(y)-5),(length(y)-5));


% background = (xmax-xmin)*beta(5) + 1/2.*beta(4).*(xmax^2-xmin^2)
title([titleName 'med p-value: ' num2str(pValue)])
counts = abs(2*sqrt(pi)*beta(2)*beta(3));
countUns = sqrt(counts+(2*sqrt(pi)*beta(3)).^2*us(2,2).^2+(2*sqrt(pi)*beta(2)).^2*us(3,3).^2);
end
function y = fitfunction(beta,x)
    y = beta(6).*x.^2+beta(4).*x+beta(5)+abs(beta(3)).*exp(-((x-beta(1))./(2.*abs(beta(2)))).^2);
end