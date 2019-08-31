clear all; close all;clc;
filename = 'F:\skole\UNI\5. Semester\Experiment 3\ExperimentalPhysics3\Kalibrering_Cs137_ch001.txt';
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

%%
% x = linspace(0,100,100);
% y = 500.*normpdf(x,50,10)+normrnd(0,1,size(x));
% yerr = ones(size(x));
figure
errorbar(X,Y,Yerr,'.')
hold on
xlabel('Channel')
ylabel('Counts')

name = input('Name of spectrum:')
title(name)
n = input('Number of peaks you want to input:')
for i = 1:n
    value(i) =input('Value of peak:')

    [x2,y2] = ginput(1);
    plot(x2,y2,'*')
    [x3,y3] = ginput(1);
    plot(x3,y3,'*')

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

beta0 = [(x2+x3)/2,(x3-x2)/3,max(y),20];
plot(x,fitfunction(beta0,x))
w = 1./yerr.^2;
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x,y,@fitfunction,beta0,'weights',w);
plot(x,fitfunction(beta,x))
us = CovB/MSE;
MSE
pValue(i) = 1-chi2cdf(MSE*(length(y)-4),(length(y)-4));


plot([beta(1)+us(1,1),beta(1)+us(1,1)],[0,max(y)])
plot([beta(1)-us(1,1),beta(1)-us(1,1)],[0,max(y)])

peakChannel(i) = beta(1,1);
peakUns(i) = us(1,1);


end

function y = fitfunction(beta,x)
    y = beta(4)+beta(3).*exp(-((x-beta(1))./(2.*beta(2))).^2);
end