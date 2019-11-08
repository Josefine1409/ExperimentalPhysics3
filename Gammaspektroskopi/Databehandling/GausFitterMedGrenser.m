clear all; close all;clc;
filenames = {'Co60_ch000.txt','Cs137_ch000.txt','Ra226_ch000.txt'};
names = {'Co60','Cs137','Ra226'};

peakValues = {[1.1732,1.3325],...
    [0.6617],...
    [0.35193,...
    0.60931,...
    0.76836,...
    1.1203,....
    1.7645]};
peakBorders = {[905,1019;1019,1151],....
    [500,600],...
    [259,326;...
    436,533;...
    551,657;...
    809,950;...
    1293,1461]};


Es = linspace(0.1,1.5,1000);
linFun =@(beta,x) (x-beta(2))/beta(1);

data = [];
for i = 1:length(filenames)
    [X,Y,Yerr] = hisFraData(filenames{i});
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


function [X,Y,Yerr] = hisFraData(filename)
filename
addpath('..\Energy calibration\')
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