clear all; close all;clc;
linFun =@(beta,x) (x-beta(2))./beta(1);

theta = 160/180*pi;
E  = [350,400]
HMasses = [1.0078,2.0141,3.0160]
K2 = @(m1,m2) ((m1.*cos(theta)+sqrt(m2.^2-m1.^2*sin(theta)^2))./(m1+m2)).^2
KC  = K2(HMasses,ones(1,3).*12.0107 )./[1,2,3]
KAu = K2(HMasses,ones(1,3).*196.96657 )./[1,2,3]

% S1Ag = 
% 
% E1Ag  = 
% peakValues ={[KAu(1)*E(1),KC(1)*E(1)],[KAu(2)*E(1),KC(2)*E(1)],[KAu(3)*E(1)]...
%             ,[KAu(1)*E(2),KC(1)*E(2)],[KAu(2)*E(2),KC(2)*E(2)],[KAu(3)*E(2)]}
% 
% peakBorders = {[400,460;280,360],[180,220;122,162],[90,150],...
%     [465,520;340,380],[205,254;148,189],[117,167]};
% 
% peakValues ={[KAu(1)*E(1),KC(1)*E(1)],[KAu(2)*E(1)],[KAu(3)*E(1)]...
%             ,[KAu(1)*E(2),KC(1)*E(2)],[KAu(2)*E(2)],[KAu(3)*E(2)]}
% 
% peakBorders = {[400,460;280,360],[180,220],[90,150],...
%     [465,520;340,380],[205,254],[117,167]};



peakValues ={[KAu(1)*E(1)],[KAu(2)*E(1)],[KAu(3)*E(1)],[KAu(1)*E(2)],[KAu(2)*E(2)],[KAu(3)*E(2)]}
peakBorders = {[400,460],[180,220],[90,150],...
    [465,520],[205,254],[117,167]};

%%
filenames = {};
names = {};
Es = {'350 KeV','400 KeV'}
data = [];

for j = 1:2
    for i = 1:3
        energy = Es{j}
%       filenames{end+1} = ['\Initialenergi ',energy,'\H',num2str(i),'+.asc'];
%         path(['..\Data\Kalibrering\Initialenergi ' energy])
%         filenames{end+1} = ['H',num2str(i),'+.asc'];
%         names{end+1} = [energy,', H',num2str(i),'+'];
        
        filename =  fullfile(['..\Data\Kalibrering\Initialenergi ' energy],['H',num2str(i),'+.asc'])
        [X,Y,Yerr] =hisFraData(filename);
        data = [data, fitGaussInSpectrum(X,Y,Yerr,[energy,', H',num2str(i),'+'],peakValues{i+(j-1)*3},peakBorders{i+(j-1)*3})];
        
        
    end
end

% 
% for i = 1:length(filenames)
%     [X,Y,Yerr] =hisFraData(filenames{i});
%     data = [data, fitGaussInSpectrum(X,Y,Yerr,names{i},peakValues{i},peakBorders{i})];
% end
Es = linspace(50,400,1000);

x =data(3,:);
y = data(1,:);
yerr = data(2,:);

beta0 = [1,2];
linFun =@(beta,x) (x-beta(2))./beta(1);
w = 1./yerr.^2;
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x,y,linFun,beta0,'weights',w);

[Ypred,delta] = nlpredci(linFun,Es,beta,R,'jacobian',J,'alpha',0.35);
ci = nlparci(beta,R,'jacobian',J,'alpha',0.35)

figure
hold on
xlabel('Energy (E) [keV]')
ylabel('Channel number (Ch)')
set(gca,'FontSize',15) 

plot(Es,Ypred,'-','linewidth',1)
errorbar(x,y',yerr,'.','markersize',8)
legend('Fit','Data','Location','northwest')


figure
xlabel('Energy (E) [keV]')
ylabel('Residual Channel number (\DeltaCh)')
set(gca,'FontSize',15) 
hold on

plot(Es,0.*Es,'b','linewidth',1)

[YpredMeasurement,deltaM] = nlpredci(@(beta,x) linFun(beta,x),x,beta,R,'jacobian',J,'alpha',0.35);
errorbar(x,y-YpredMeasurement',yerr,'.r','markersize',8)
plot(Es,delta,'k--','linewidth',1)
plot(Es,-delta,'k--','linewidth',1)

legend('Fit','Data','Fit confidence','Location','southwest')

us = CovB/MSE
MSE
pValue = 1-chi2cdf(MSE*(length(y)-2),(length(y)-2))

disp(['Energi som funktion af channel: E= ' num2str(beta(1)) '+-' num2str((ci(1,2)-ci(1,1))/2) 'Channel+' num2str(beta(2)) '+-' num2str((ci(2,2)-ci(2,1))/2) 'MeV'])


%%
function [X,Y,Yerr] = hisFraData(filename)
delimiter = ' ';
startRow = 12;

formatSpec = '%f%*s%*s%[^\n\r]';
fileID = fopen(filename);
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
Spectrum = dataArray{:, 1};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

X = 1:length(Spectrum);
Y = Spectrum;
Yerr = sqrt(Y) +(Y==0);

end

%%
% data har sturktur [peakChannel,peakUns,peakValue]

function data = fitGaussInSpectrum(X,Y,Yerr,name,peakValue,peakBorder)
n = length(peakValue);

figure
errorbar(X*0.76468+12.393,Y,Yerr,'.')
xlabel('Energy (E) [keV]')
ylabel('Counts (n)')
set(gca,'FontSize',15) 
xlim([0,400])

hold on
for i = 1:n
    plot([peakValue(i),peakValue(i)],[0,500])
end
figure
histogram('BinEdges',[0,X(1:end)],'BinCounts',Y,'EdgeColor','none')
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
w = 1./yerr.^2;
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x',y,@fitfunction,beta0,'weights',w);
beta(1);
beta(2);
plot(x,fitfunction(beta,x),'k--','linewidth',2)
xlim([0,500])
us = CovB/MSE;
mse =MSE;
MSECount(i) = MSE;
pValue(i) = 1-chi2cdf(MSE*(length(y)-5),(length(y)-5));
P_Value = 1-chi2cdf(MSE*(length(y)-5),(length(y)-5));


txt = text(beta(1),max(y)+10,['\leftarrow' num2str(peakValue(i)) ' keV']);
set(txt,'Rotation',90);
set(txt,'FontSize',14);

peakChannel(i) = beta(1,1);
peakUns(i) = us(1,1);
ci = nlparci(beta,R,'jacobian',J,'alpha',0.35);
peakUns(i) = (ci(1,2)-ci(1,1))/2;

end

peakUns;
peakChannel;
peakValue;
data = [peakChannel;peakUns;peakValue;pValue;MSECount];
end


function y = fitfunction(beta,x)
    y = beta(4).*x+beta(5)+beta(3).*exp(-((x-beta(1))./(beta(2))).^2./2);
end