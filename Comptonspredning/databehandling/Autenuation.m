clear all; close all;clc;
folderName = {'Aluminum'}
name = {'Al'}
channelInterval= [[520,580]]
plateThikness = {[0,10.11,20.22,30.34,40.46,50.55,60.65,70.79,80.9]}

for i = 1:length(name)
    pathStart = ['..\data\AttenuationCoefficient\' folderName{i} '\AttenuationCoefficient_' name{i} '_'];
        for j = 1:length(plateThikness{i})
            path  = [pathStart num2str(j-1) 'Plates_ch001.txt']
            [X,Y,Yerr] = hisFraData(path);
            counts(j) = sum(Y((channelInterval(1,1):channelInterval(1,2))))
            figure
            plot(X,Y,'.')
            hold on 
            plot(X([channelInterval(1,1),channelInterval(1,2)]),Y([channelInterval(1,1),channelInterval(1,2)]),'*')
        end 
    figure
    cErr = sqrt(counts)
    x = plateThikness{i}
    y = counts
    yerr = cErr
    errorbar(x,y,yerr,'.')
    hold on
    beta0 = [y(1)-y(end),-0.1,y(end),0];
    linFun =@(beta,x) beta(1).*exp(beta(2).*x)+beta(3)+beta(4).*x
    plot(x,linFun(beta0,x))
    w = 1./yerr.^2;
    [beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x,y,@(beta,x) linFun(beta,x),beta0,'weights',w)
    hold on
    plot(x,linFun(beta,x))
    us = CovB/MSE;
    MSE
    pValue = 1-chi2cdf(MSE*(length(y)-2),(length(y)-2))
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
