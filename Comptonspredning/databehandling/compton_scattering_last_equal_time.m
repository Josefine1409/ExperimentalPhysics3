clc; clear all; close all;
%% Import Data

%Vinkel for denne m�ling.
theta=70;
n=1;
% path1 = ['C:\Users\Lukas\CloudStation\skole\UNI\5. Semester\Experiment 3\ExperimentalPhysics3\Comptonspredning\data\GivenData\',num2str(theta(n)),'deg-30min_ch000.txt'];

path1 = 'C:\Users\Lukas\CloudStation\skole\UNI\5. Semester\Experiment 3\ExperimentalPhysics3\Comptonspredning\data\GivenData\70deg-30min_ch000.txt';


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
path2 = 'C:\Users\Lukas\CloudStation\skole\UNI\5. Semester\Experiment 3\ExperimentalPhysics3\Comptonspredning\data\GivenData\70deg-30min_ch001.txt';
delimiter = ' ';
startRow = 6;
formatSpec = '%f%f%f%[^\n\r]';
fileID = fopen(path2,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
timestamp2 = dataArray{:, 1};
channel2 = dataArray{:, 2};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;




%M�linger. M1 er f�rste m�ling af spredningen, M2 er anden m�ling af
%spredningen. F�rste s�jle er tid, anden s�jle er energi.
M1=[timestamp1 channel1];

M2=[timestamp2 channel2];

% Nr=1e6
% M1=M1(1:Nr,:);
% 
% M2=M2(1:Nr,:);

%M�ledata M1's tider fittes line�rt til indekser, for at begr�nse data
%indenfor tidstolerancen. Indekstolerancen s�rger for, at vi ikke misser
%data udenfor det pr�cise, line�re fit grundet spredning.

LM1=length(M1(:,1));
 a=nlinfit([1:length(M1(:,1))]',M1(:,1),@(a,x)a(1)*x+a(2),[1 0]);

%% Tolerance
close all
LM2=length(M2);

%Tolerancer og gamma-energier fra kilden. Tidstolerancen m�ler hvor lang
%tid bagud, vi vil acceptere en m�ling. Energitolerancen beskriver, hvor
%meget energibevarelsen m� afvige enten positivt eller negativt.
% TolTid=0.25/3;
TolTid=70;

%S�rger for at man ikke misser gode m�linger udenfor intervallet grundet
%spredning fra det line�re fit.
IndeksTol=1000;


data=[];
for i=1:length(M2(:,1))
%     Perc=100*i/length(M2(:,1))
    
    %Indeksinterval, hvor tidskravet er opfyldt.
%     M2(i,:)
    time1_less_than_time2 = 0;
    currentGuess = round((M2(i,1)-TolTid-a(2))/a(1));
    while time1_less_than_time2 
        currentGuess  = currentGuess-1;
        time1_less_than_time2 = M1(currentGuess,1) <M2(i,1);
    end

    I1=max(round((M2(i,1)-TolTid-a(2))/a(1))-IndeksTol,1);
    I2=min(floor((M2(i,1)-a(2))/a(1))+IndeksTol,LM1);

    %Begr�nsning af m�linger til det relevante interval.
%     M1=M1(I1-sumI1:length(M1),[1 2]);
%     m1=M1(1:I2-sumI1,[1 2]);
% 
%     sumI1=I1-1;

m1=M1(I1:I2,[1 2]);


    %Pr�csise krav ift. tidstolerancen overholdes.

    m1=m1(m1(:,1)<=M2(i,1),[1 2]);
    m1=m1(m1(:,1)>=M2(i,1)-TolTid,[1 2]);
    
    m1(:,1)=M2(i,1)-m1(:,1);
    
   data=[data;m1];
end

histogram(data(:,1),20)
figure
histogram(data(:,2),20,'FaceColor','r')
figure
histogram(M2(:,2),20,'FaceColor','b','binlimits',[0 900])
figure
histogram([M2(:,2);data(:,2)],20,'binlimits',[0 900])

%% Sort Data
tic

%Tolerancer og gamma-energier fra kilden. Tidstolerancen m�ler hvor lang
%tid bagud, vi vil acceptere en m�ling. Energitolerancen beskriver, hvor
%meget energibevarelsen m� afvige enten positivt eller negativt.
% TolTid=0.25/3;
TolTid=20;
TolEnergi=200;
E_Cs=548.6;

%S�rger for at man ikke misser gode m�linger udenfor intervallet grundet
%spredning fra det line�re fit.
IndeksTol=100;

N=0;
% plot(M1(:,1),'r.')
% hold on
% plot([1:length(M1(:,1))]',a(1)*[1:length(M1(:,1))]'+a(2))
sumI1=0;

LM1=length(M1(:,1));

for i=1:length(M2(:,1))
%     Perc=100*i/length(M2(:,1))
    
    %Indeksinterval, hvor tidskravet er opfyldt.
%     M2(i,:)

    I1=max(round((M2(i,1)-TolTid-a(2))/a(1))-IndeksTol,1);
    I2=min(floor((M2(i,1)-a(2))/a(1))+IndeksTol,LM1);

    %Begr�nsning af m�linger til det relevante interval.
%     M1=M1(I1-sumI1:length(M1),[1 2]);
%     m1=M1(1:I2-sumI1,[1 2]);
% 
%     sumI1=I1-1;

m1=M1(I1:I2,[1 2]);


    %Pr�csise krav ift. tidstolerancen overholdes.

    m1=m1(m1(:,1)<=M2(i,1),[1 2]);
    m1=m1(m1(:,1)>=M2(i,1)-TolTid,[1 2]);
    
    
    %Energitolerancen checkes indenfor relevant data.
   
    m1=m1(m1(:,2)<E_Cs-M2(i,2)+TolEnergi & m1(:,2)>E_Cs-M2(i,2)-TolEnergi,[1 2]);
    
    
    %Hvis der er gyldige m�linger, t�ller vi, hvis ikke, t�ller vi ikke.
    
    N=N+(length(m1(:,1))>0);
    
end

Antal_godkendte_maalinger=N
Vinkel=theta(n)   
toc