clc; clear all; close all;

%Vinkel for denne m�ling.
theta=;


%M�linger. M1 er f�rste m�ling af spredningen, M2 er anden m�ling af
%spredningen. F�rste s�jle er tid, anden s�jle er energi.
M1=[];

M2=[];

%Tolerancer og gamma-energier fra kilden. Tidstolerancen m�ler hvor lang
%tid bagud, vi vil acceptere en m�ling. Energitolerancen beskriver, hvor
%meget energibevarelsen m� afvige enten positivt eller negativt.
TolTid=0.3;
TolEnergi=1;
E_Cs=0.5;

%S�rger for at man ikke misser gode m�linger udenfor intervallet grundet
%spredning fra det line�re fit.
IndeksTol=10;

%M�ledata M1's tider fittes line�rt til indekser, for at begr�nse data
%indenfor tidstolerancen. Indekstolerancen s�rger for, at vi ikke misser
%data udenfor det pr�cise, line�re fit grundet spredning.
a=nlinfit([1:length(M1(:,1))]',M1(:,1),@(a,x)a(1)*x+a(2),[1 0]);
N=0;

for i=1:length(M2(:,1))
    
    %Indeksinterval, hvor tidskravet er opfyldt.
    I1=max(round((M2(i,1)-TolTid-a(2))/a(1))+IndeksTol,1);
    I2=min(floor((M2(i,1)-a(2))/a(1))-IndeksTol,length(M1(:,1)));
    
    %Begr�nsning af m�linger til det relevante interval.
    m=M1(I1:I2,[1 2]);
    
    %Pr�csise krav ift. tidstolerancen overholdes.
    m=m(m(:,1)<M2(i,1),[1 2]);
    m=m(m(:,1)>M2(i,1)-TolTid,[1 2]);
    
    %Energitolerancen checkes indenfor relevant data.
    m=m(m(:,2)<E_Cs+TolEnergi & m(:,2)>E_Cs-TolEnergi,[1 2]);
    
    %Hvis der er gyldige m�linger, t�ller vi, hvis ikke, t�ller vi ikke.
    N=N+(length(m(:,1))>0);
end

[Antal_godkendte_maalinger Vinkel]=[N theta]