clc; clear all; close all;

%Vinkel for denne måling.
theta=20:10:120;
n=5;
path1 = ['C:\Users\Rasmu\Documents\Eksperimentel 3\Git\ExperimentalPhysics3\Comptonspredning\data\GivenData\',num2str(theta(n)),'deg-30min_ch000.txt'];

delimiter = ' ';
startRow = 6;
formatSpec = '%f%f%f%[^\n\r]';
fileID = fopen(path1,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
timestamp1 = dataArray{:, 1};
channel1 = dataArray{:, 2};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

path2 = ['C:\Users\Rasmu\Documents\Eksperimentel 3\Git\ExperimentalPhysics3\Comptonspredning\data\GivenData\',num2str(theta(n)),'deg-30min_ch001.txt'];

delimiter = ' ';
startRow = 6;
formatSpec = '%f%f%f%[^\n\r]';
fileID = fopen(path2,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
timestamp2 = dataArray{:, 1};
channel2 = dataArray{:, 2};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;



%Målinger. M1 er første måling af spredningen, M2 er anden måling af
%spredningen. Første søjle er tid, anden søjle er energi.
M1=[timestamp1 channel1];

M2=[timestamp2 channel2];

Nr=1e5;
M1=M1(1:Nr,:);

M2=M2(1:Nr,:);

%Tolerancer og gamma-energier fra kilden. Tidstolerancen måler hvor lang
%tid bagud, vi vil acceptere en måling. Energitolerancen beskriver, hvor
%meget energibevarelsen må afvige enten positivt eller negativt.
TolTid=0.25/3;
TilTid=20;
TolEnergi=14.5;
E_Cs=548.6;

%Sørger for at man ikke misser gode målinger udenfor intervallet grundet
%spredning fra det lineære fit.
IndeksTol=;

%Måledata M1's tider fittes lineært til indekser, for at begrænse data
%indenfor tidstolerancen. Indekstolerancen sørger for, at vi ikke misser
%data udenfor det præcise, lineære fit grundet spredning.
a=nlinfit([1:length(M1(:,1))]',M1(:,1),@(a,x)a(1)*x+a(2),[1 0]);
N=0;
% plot(M1(:,1),'r.')
% hold on
% plot([1:length(M1(:,1))]',a(1)*[1:length(M1(:,1))]'+a(2))

for i=1:length(M2(:,1))
%     Perc=100*i/length(M2(:,1))
    
    %Indeksinterval, hvor tidskravet er opfyldt.
%     M2(i,:)
    I1=max(round((M2(i,1)-TolTid-a(2))/a(1))-IndeksTol,1);
    I2=min(floor((M2(i,1)-a(2))/a(1))+IndeksTol,length(M1(:,1)));
    
    %Begrænsning af målinger til det relevante interval.
    m1=M1(I1:I2,[1 2]);

    %Præcsise krav ift. tidstolerancen overholdes.
    m1=m1(m1(:,1)<=M2(i,1),[1 2]);
    m1=m1(m1(:,1)>=M2(i,1)-TolTid,[1 2]);
    
    %Energitolerancen checkes indenfor relevant data.
    m1=m1(m1(:,2)<E_Cs-M2(i,2)+TolEnergi & m1(:,2)>E_Cs-M2(i,2)-TolEnergi,[1 2]);
    
    %Hvis der er gyldige målinger, tæller vi, hvis ikke, tæller vi ikke.
    N=N+(length(m1(:,1))>0);
end

Antal_godkendte_maalinger=N
Vinkel=theta(n)   