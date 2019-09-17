%% Import Data
clc; clear all; close all;

thetas=[70 110];
times=[15];

%Vinkel for denne måling (n).
n=1;

theta=thetas(n);
time=times(n)

% path1 = ['C:\Users\Rasmu\Documents\Eksperimentel 3\Git\ExperimentalPhysics3\Comptonspredning\data\GivenData\',num2str(theta(n)),'deg-30min_ch000.txt'];
path1 = ['C:\Users\Rasmu\Documents\Eksperimentel 3\Git\ExperimentalPhysics3\Comptonspredning\data\Compton_',num2str(theta),'deg_',num2str(time),'min_ch000.txt'];

delimiter = ' ';
startRow = 6;
formatSpec = '%f%f%f%[^\n\r]';
fileID = fopen(path1,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
timestamp1 = dataArray{:, 1};
channel1 = dataArray{:, 2};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

% path2 = ['C:\Users\Rasmu\Documents\Eksperimentel 3\Git\ExperimentalPhysics3\Comptonspredning\data\GivenData\',num2str(theta(n)),'deg-30min_ch001.txt'];
path2 = ['C:\Users\Rasmu\Documents\Eksperimentel 3\Git\ExperimentalPhysics3\Comptonspredning\data\Compton_70deg_15min_ch001.txt'];

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

% Nr=1e6
% M1=M1(1:Nr,:);
% 
% M2=M2(1:Nr,:);

%Måledata M1's tider fittes lineært til indekser, for at begrænse data
%indenfor tidstolerancen. Indekstolerancen sørger for, at vi ikke misser
%data udenfor det præcise, lineære fit grundet spredning.
a=nlinfit([1:length(M1(:,1))]',M1(:,1),@(a,x)a(1)*x+a(2),[1 0]);

D=zeros(length(M1(:,1)),1);
for i=1:length(M1(:,1))
    D(i)=norm(i-(M1(i,1)-a(2))/a(1));
end
%Sørger for at man ikke misser gode målinger udenfor intervallet grundet
%spredning fra det lineære fit.
IndeksTol=ceil(max(D));

%% Tolerance
close all
LM2=length(M2);
LM1=length(M1(:,1));

%Tolerancer og gamma-energier fra kilden. Tidstolerancen måler hvor lang
%tid bagud, vi vil acceptere en måling. Energitolerancen beskriver, hvor
%meget energibevarelsen må afvige enten positivt eller negativt.
% TolTid=0.25/3;
TolTid=40;

data=[];
data1=[];
m2=zeros(length(M2),1);
for i=1:length(M2(:,1))
%     Perc=100*i/length(M2(:,1))
    
    %Indeksinterval, hvor tidskravet er opfyldt.
%     M2(i,:)

    I1=max(round((M2(i,1)-TolTid-a(2))/a(1))-IndeksTol,1);
    I2=min(floor((M2(i,1)-a(2))/a(1))+IndeksTol,LM1);

    %Begrænsning af målinger til det relevante interval.
%     M1=M1(I1-sumI1:length(M1),[1 2]);
%     m1=M1(1:I2-sumI1,[1 2]);
% 
%     sumI1=I1-1;

m1=M1(I1:I2,[1 2]);


    %Præcsise krav ift. tidstolerancen overholdes.

    m1=m1(m1(:,1)<=M2(i,1),[1 2]);
    m1=m1(m1(:,1)>=M2(i,1)-TolTid,[1 2]);
    
    m1(:,1)=M2(i,1)-m1(:,1);
    data1=[data1 m1(:,2)];
    m1(:,2)=M2(i,2)+m1(:,2);
    
    m2(i)=(M2(i,2)+50)*(1-isempty(m1(:,1)));
    
   data=[data;m1];
end
m2=m2-50;
% length(data)
histogram(data(:,1))
figure
histogram(data(:,2),'binlimits',[0 1200])
figure
histogram(data1,'binlimits',[0 1200])
figure
histogram(m2,'binlimits',[0 1200])

%% Sort Data

%Tolerancer og gamma-energier fra kilden. Tidstolerancen måler hvor lang
%tid bagud, vi vil acceptere en måling. Energitolerancen beskriver, hvor
%meget energibevarelsen må afvige enten positivt eller negativt.
% TolTid=0.25/3;
TolTid=30;
TolEnergi=130;
E_Cs=870;


N=0;
% plot(M1(:,1),'r.')
% hold on
% plot([1:length(M1(:,1))]',a(1)*[1:length(M1(:,1))]'+a(2))
sumI1=0;



for i=1:length(M2(:,1))
%     Perc=100*i/length(M2(:,1))
    
    %Indeksinterval, hvor tidskravet er opfyldt.
%     M2(i,:)

    I1=max(round((M2(i,1)-TolTid-a(2))/a(1))-IndeksTol,1);
    I2=min(floor((M2(i,1)-a(2))/a(1))+IndeksTol,LM1);

    %Begrænsning af målinger til det relevante interval.
%     M1=M1(I1-sumI1:length(M1),[1 2]);
%     m1=M1(1:I2-sumI1,[1 2]);
% 
%     sumI1=I1-1;

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
Vinkel=theta
