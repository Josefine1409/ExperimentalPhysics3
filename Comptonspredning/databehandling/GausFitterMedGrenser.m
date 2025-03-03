clear all; close all;clc;
% filenames = {'Kalibrering_Co60_ch001.txt','Kalibrering_Cs137_ch001.txt','Kalibrering_Ra226_ch001.txt'};
% names = {'Co60','Cs137','Ra226'};
% 
% peakValues = {[1.1732,1.3325],[0.6617],[0.35193,0.60931,0.76836,1.1203,1.7645]};
% peakBorders = {[905,1019;1019,1151],[500,600],[259,326;436,533;551,657;809,950;1293,1461]};
% 

% 
% filenames = {'Kalibrering_Ra226_300s_ch001.txt'};
% names = {'Ra226'};
% 
% peakValues = {[0.1861,0.24198,0.29522,0.35193,0.60931,1.1203 ]};
% peakBorders = {[137,168;183,222;222,263;263,316;440,540;834,942]};
% 
% linFun =@(beta,x) (x-beta(2))/beta(1);
% 
% data = [];
% for i = 1:length(filenames)
%     [X,Y,Yerr] =hisFraData(filenames{i});
%     data = [data, fitGaussInSpectrum(X,Y,Yerr,names{i},peakValues{i},peakBorders{i})];
% end
% x =data(3,:);
% y = data(1,:);
% yerr = data(2,:);
% 
% beta0 = [1,2];
% linFun =@(beta,x) (x-beta(2))/beta(1);
% w = 1./yerr.^2;
% [BETA,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x,y,@(beta,x) linFun(beta,x),beta0,'weights',w);
% 
% 
% 
% 
% 
% 
% 
Es = linspace(0.1,1.5,1000);
% 
% 
% %%
% filenames = {'Kalibrering_Co60_ch001.txt','Kalibrering_Cs137_ch001.txt'};
% names = {'Co60','Cs137'};
% 
% peakValues = {[1.1732,1.3325],[0.66164]};
% peakBorders = {[905,1019;1019,1151],[500,600]};
% 
% Data = [];
% 
% for i = 1:length(filenames)
%     [X,Y,Yerr] =hisFraData(filenames{i});
%     Data = [Data, fitGaussInSpectrum(X,Y,Yerr,names{i},peakValues{i},peakBorders{i})];
% end
% figure(10)
% hold on
% x =Data(3,:);
% y = Data(1,:);
% yerr = Data(2,:);
% 
% beta0 = [0.0012668,-0.013154];
% 
% w = 1./yerr.^2;
% [beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x,y,@(beta,x) linFun(beta,x),beta0,'weights',w);
% plot(Es,linFun(beta,Es)-linFun(BETA,Es),'r','linewidth',1)
% % plot(Es,linFun(beta0,Es),'r','linewidth',1)
% 
% hold on
% 
% errorbar(x,y-linFun(BETA,x),yerr,'.b','markersize',8);
% 
% us = CovB/MSE
% mse = MSE
% pValue = 1-chi2cdf(MSE*(length(y)-2),(length(y)-2))
% 
% disp(['Energi som funktion af channel: E= ' num2str(beta(1)) '+-' num2str(us(1,1)) 'Channel-' num2str(beta(2)) '+-' num2str(us(2,2)) 'MeV'])
% % xlim([data(3,1)-0.05,data(3,end)+0.05])
% 
% 
% 
% %%
% filenames = {'Kalibrering_Ra226_ch001.txt'};
% names = {'Ra226'};
% 
% peakValues = {[0.35193,0.60931,0.76836,1.1203,1.3777 ,1.7645]};
% peakBorders = {[259,326;436,533;551,657;809,950;1040,1120;1293,1461]};
% 
% [X,Y,Yerr] = hisFraData(filenames{1});
% data = fitGaussInSpectrum(X,Y,Yerr,names{1},peakValues{1},peakBorders{1});
% figure(10)
% hold on
% xlabel('Energy [MeV]')
% ylabel('Channel number')
% set(gca,'FontSize',15) 
% 
% x =data(3,:);
% y = data(1,:);
% yerr = data(2,:);
% 
% w = 1./yerr.^2;
% [beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x,y,@(beta,x) linFun(beta,x),beta0,'weights',w);
% hold on
% plot(Es,linFun(beta,Es)-linFun(BETA,Es),'b','linewidth',1)
% errorbar(x,y-linFun(BETA,x),yerr,'.r','markersize',8)
% 
% us = CovB/MSE
% MSE
% pValue = 1-chi2cdf(MSE*(length(y)-2),(length(y)-2))
% 
% disp(['Energi som funktion af channel: E= ' num2str(beta(1)) '+-' num2str(us(1,1)) 'Channel-' num2str(beta(2)) '+-' num2str(us(2,2)) 'MeV'])



%%

    
filenames = {'Kalibrering_Ra226_300s_ch001.txt'};
% peakValues = {[0.1861,0.24198,0.29521,0.35192,0.60931]};
peakValues = {[0.1861,0.24198,0.29522,0.35193,0.60931,1.1203 ]};
 names = {'Ra226'};
peakBorders = {[137,168;183,222;222,263;263,316;440,540;834,942]};
% peakValues = {[0.1861,0.24198,0.29522,0.35193,0.60931]};
% 
% peakBorders = {[137,168;183,222;222,263;263,316;440,540]};

linFun =@(beta,x) (x-beta(2))/beta(1);
data = [];

for i = 1:length(filenames)
    [X,Y,Yerr] =hisFraData(filenames{i});
    data = [data, fitGaussInSpectrum(X,Y,Yerr,names{i},peakValues{i},peakBorders{i})];
end



x =data(3,:);
y = data(1,:);
yerr = data(2,:);

beta0 = [1,2];
linFun =@(beta,x) (x-beta(2))/beta(1);
w = 1./yerr.^2;
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x,y,@(beta,x) linFun(beta,x),beta0,'weights',w);

[Ypred,delta] = nlpredci(@(beta,x) linFun(beta,x),Es,beta,R,'jacobian',J,'alpha',0.35);
ci = nlparci(beta,R,'jacobian',J,'alpha',0.35)

figure(10)
hold on
xlabel('Energy (E) [MeV]')
ylabel('Channel number (Ch)')
set(gca,'FontSize',15) 

plot(Es,Ypred,'-','linewidth',1)
errorbar(x,y',yerr,'.','markersize',8)

% plot(Es,Ypred+delta,'b--','linewidth',1)
% plot(Es,Ypred-delta,'b--','linewidth',1)

legend('Fit','Data','Location','northwest')



figure(11)
xlabel('Energy (E) [MeV]')
ylabel('Residual Channel number (\DeltaCh)')
set(gca,'FontSize',15) 
hold on

plot(Es,0.*Es,'b','linewidth',1)

[YpredMeasurement,deltaM] = nlpredci(@(beta,x) linFun(beta,x),x,beta,R,'jacobian',J,'alpha',0.35);
errorbar(x,y-YpredMeasurement',yerr,'.r','markersize',8)
plot(Es,delta,'k--','linewidth',1)
plot(Es,-delta,'k--','linewidth',1)

legend('Fit','Data','Fit confidence','Location','southwest')




% plot(Es,linFun(beta,Es),'g','linewidth',1)
% plot(Es,linFun(beta,Es)-linFun(BETA,Es)+delta),'g--','linewidth',1)
% plot(Es,linFun(beta,Es)-linFun(BETA,Es)-delta),'g--','linewidth',1)

% plot(Es,Ypred,'g','linewidth',1)
% plot(Es,Ypred+delta,'g--','linewidth',1)
% plot(Es,Ypred-delta,'g--','linewidth',1)

% plot(Es,0.*Es,'g','linewidth',1)
% plot(Es,delta,'g--','linewidth',1)
% plot(Es,-delta,'g--','linewidth',1)
% 
% 
% % errorbar(x,y-linFun(BETA,x),yerr,'.k','markersize',8)
% [YpredMeasurement,deltaM] = nlpredci(@(beta,x) linFun(beta,x),x,beta,R,'jacobian',J,'alpha',0.35);
% errorbar(x,y-YpredMeasurement',yerr,'.k','markersize',8)

us = CovB/MSE
MSE
pValue = 1-chi2cdf(MSE*(length(y)-2),(length(y)-2))

disp(['Energi som funktion af channel: E= ' num2str(beta(1)) '+-' num2str((ci(1,2)-ci(1,1))/2) 'Channel+' num2str(beta(2)) '+-' num2str((ci(2,2)-ci(2,1))/2) 'MeV'])


% legend('Fit 1. day','Data 1. day','Fit 2. day','Data 2. day','Fit 3. day','Data 3. day','Location','northwest')

%%


function [X,Y,Yerr] = hisFraData(filename)
addpath('..\data\Kalibrering')
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
histogram
histcounts
%%
% data har sturktur [peakChannel,peakUns,peakValue]

function data = fitGaussInSpectrum(X,Y,Yerr,name,peakValue,peakBorder)
n = length(peakValue);


figure
errorbar(X*0.0012668-0.013154,Y,Yerr,'.')
% errorbar(X,Y,Yerr,'.','markersize',10)

xlabel('Energy (E) [MeV]')
% xlabel('Channel number')
ylabel('Counts (n)')
set(gca,'FontSize',15) 

hold on
for i = 1:n
    plot([peakValue(i),peakValue(i)],[0,500])
end
figure
% histogram('Categories',X,'BinCounts',Y)
histogram('BinEdges',[1/2,X(1:end)+1/2],'BinCounts',Y,'EdgeColor','none')
% errorbar(X,Y,Yerr,'.','markersize',10)
hold on
xlabel('Channel number (Ch)')
ylabel('Counts (n)')
set(gca,'FontSize',15) 
title(name)

for i = 1:n

x2 =   peakBorder(i,1);
x3 =   peakBorder(i,2);  

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

beta0 = [(x2+x3)/2,(x3-x2)/3,max(y),0,20];
% plot(x,fitfunction(beta0,x))
w = 1./yerr.^2;
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x,y,@fitfunction,beta0,'weights',w);
beta(1);
beta(2);
plot(x,fitfunction(beta,x),'k--','linewidth',2)
% plot(x,beta(4).*x+beta(5))
us = CovB/MSE;
mse =MSE;
MSECount(i) = MSE;
pValue(i) = 1-chi2cdf(MSE*(length(y)-5),(length(y)-5));
P_Value = 1-chi2cdf(MSE*(length(y)-5),(length(y)-5));


txt = text(beta(1),max(y)+10,['\leftarrow' num2str(peakValue(i)) ' MeV']);
set(txt,'Rotation',90);
set(txt,'FontSize',14);

% plot([beta(1)+us(1,1),beta(1)+us(1,1)],[0,max(y)])
% plot([beta(1)-us(1,1),beta(1)-us(1,1)],[0,max(y)])

peakChannel(i) = beta(1,1);
peakUns(i) = us(1,1);
ci = nlparci(beta,R,'jacobian',J,'alpha',0.35)
peakUns(i) = (ci(1,2)-ci(1,1))/2;
% peakUns(i) = us(1,1);

end

peakUns;
peakChannel;
peakValue;
data = [peakChannel;peakUns;peakValue;pValue;MSECount];
end


function y = fitfunction(beta,x)
    y = beta(4).*x+beta(5)+beta(3).*exp(-((x-beta(1))./(beta(2))).^2./2);
end