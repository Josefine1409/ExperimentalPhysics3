clear all; close all;clc;
filenames = {'Kalibrering_Co60_ch001.txt','Kalibrering_Cs137_ch001.txt'};
names = {'Co60','Cs137'};

peakValues = {[1.1732,1.3325],[0.6617]}

data = [];
for i = 1:length(filenames)
    [X,Y,Yerr] =hisFraData(filenames{i})
    data = [data, fitGaussInSpectrum(X,Y,Yerr,names{i},peakValues{i})]
end
figure
hold on
xlabel('Energy')
ylabel('Channel')
x =data(3,:)
y = data(1,:)
yerr = data(2,:)

errorbar(x,y,yerr,'.')
beta0 = [1,2];
linFun =@(beta,x) beta(1).*x+beta(2)
w = 1./yerr.^2;
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x,y,@(beta,x) linFun(beta,x),beta0,'weights',w)
hold on
plot(x,linFun(beta,x))
us = CovB/MSE
MSE
pValue = 1-chi2cdf(MSE*(length(y)-2),(length(y)-2))

disp(['Energi som funktion af channel: E= ' num2str(beta(1)^-1) 'Channel-' num2str(beta(2)/beta(1)) 'MeV'])

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

length(channel)
X = 1:max(channel);
for i = X
    Y(i) = sum(channel==i);
end
Yerr = sqrt(Y) +(Y==0);

end

%%
% data har sturktur [peakChannel,peakUns,peakValue]

function data = fitGaussInSpectrum(X,Y,Yerr,name,peakValue)

figure
errorbar(X,Y,Yerr,'.')
hold on
xlabel('Channel')
ylabel('Counts')
title(name)
n = length(peakValue) 
for i = 1:n
    input(['Zoom in on peak:',num2str(peakValue(i))])
    
    areaChosen = 0;
    while areaChosen ~=1
        [x2,y2] = ginput(1);
        plot(x2,y2,'*')
        [x3,y3] = ginput(1);
        plot(x3,y3,'*')
        areaChosen = input('Did you choese the correct area? 1 for yes:');
    end
    
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
plot(x,fitfunction(beta0,x))
w = 1./yerr.^2;
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x,y,@fitfunction,beta0,'weights',w);
beta(1)
beta(2)
plot(x,fitfunction(beta,x))
plot(x,beta(4).*x+beta(5))
us = CovB/MSE;
MSE
MSECount(i) = MSE;
pValue(i) = 1-chi2cdf(MSE*(length(y)-5),(length(y)-5));


plot([beta(1)+us(1,1),beta(1)+us(1,1)],[0,max(y)])
plot([beta(1)-us(1,1),beta(1)-us(1,1)],[0,max(y)])

peakChannel(i) = beta(1,1);
peakUns(i) = us(1,1);
end

peakUns
peakChannel
peakValue
data = [peakChannel;peakUns;peakValue;pValue;MSECount];
end


function y = fitfunction(beta,x)
    y = beta(4).*x+beta(5)+beta(3).*exp(-((x-beta(1))./(beta(2))).^2./2);
end