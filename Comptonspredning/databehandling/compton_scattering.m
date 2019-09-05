clc; clear all; close all;

%Vinkel for denne måling.
theta=;


%Målinger. M1 er første måling af spredningen, M2 er anden måling af
%spredningen. Første søjle er tid, anden søjle er energi.
M1=[];

M2=[];

%Tolerancer og gamma-energier fra kilden. Tidstolerancen måler hvor lang
%tid bagud, vi vil acceptere en måling. Energitolerancen beskriver, hvor
%meget energibevarelsen må afvige enten positivt eller negativt.
TolTid=0.3;
TolEnergi=1;
E_Cs=0.5;

%Sørger for at man ikke misser gode målinger udenfor intervallet grundet
%spredning fra det lineære fit.
IndeksTol=10;

%Måledata M1's tider fittes lineært til indekser, for at begrænse data
%indenfor tidstolerancen. Indekstolerancen sørger for, at vi ikke misser
%data udenfor det præcise, lineære fit grundet spredning.
a=nlinfit([1:length(M1(:,1))]',M1(:,1),@(a,x)a(1)*x+a(2),[1 0]);
N=0;

for i=1:length(M2(:,1))
    
    %Indeksinterval, hvor tidskravet er opfyldt.
    I1=max(round((M2(i,1)-TolTid-a(2))/a(1))+IndeksTol,1);
    I2=min(floor((M2(i,1)-a(2))/a(1))-IndeksTol,length(M1(:,1)));
    
    %Begrænsning af målinger til det relevante interval.
    m=M1(I1:I2,[1 2]);
    
    %Præcsise krav ift. tidstolerancen overholdes.
    m=m(m(:,1)<M2(i,1),[1 2]);
    m=m(m(:,1)>M2(i,1)-TolTid,[1 2]);
    
    %Energitolerancen checkes indenfor relevant data.
    m=m(m(:,2)<E_Cs+TolEnergi & m(:,2)>E_Cs-TolEnergi,[1 2]);
    
    %Hvis der er gyldige målinger, tæller vi, hvis ikke, tæller vi ikke.
    N=N+(length(m(:,1))>0);
end

[Antal_godkendte_maalinger Vinkel]=[N theta]