clear all; close all;clc;

c2E = @(x) (x*0.76535+12.393)./1000;
E0 = 0.349;


mAu=196.966-4.4858e-4*79-0.03343120468;
mC=12-4.4858e-4*12;
mH=1.007276;
theta = 160/180*pi;
K2=@(theta,m2)((mH*cos(theta)+sqrt(m2^2-mH^2*sin(theta).^2))./(mH+m2)).^2;


I = importfile('..\Data\Thickness\Sample1_0degree_480sec.asc');
E = 1:length(I);
Ierr = sqrt(I);
figure
peakBorders = [242,353;353,437];
name = 'samp1 pos 2';
data = fitGaussInSpectrum(E,I,Ierr,name,2,peakBorders);
E2C = c2E(data(1,1));
E2CErr = (data(2,1))/0.76535/1000;
E2Au = c2E(data(1,2));
E2AuErr = (data(2,2))/0.76535/1000;

I = importfile('..\Data\Thickness\Sample1_180degree_480sec.asc');
E = 1:length(I);
Ierr = sqrt(I);

peakBorders = [242,353;380,475];
name = 'samp1 pos 1';
data = fitGaussInSpectrum(E,I,Ierr,name,2,peakBorders);
legend('show')
E1C = c2E(data(1,1));
E1CErr = (data(2,1))/0.76535/1000;
E1Au = c2E(data(1,2));
E1AuErr = (data(2,2))/0.76535/1000;
%%
SAu = @(E) stoppingpowerAu(E);
SC = @(E) stoppingpowerC(E);
a = SC(E0)*K2(theta,mAu)+SC(E2Au)/cos(pi-theta);
b = 1/(2*cos(pi-theta))*(SAu(E2Au)-SAu(E1Au));
c = 1/(2*cos(pi-theta))*(SC(E1C)-SC(E2C));
d = SAu(E0)*K2(theta,mC)+SAu(E1C)/cos(pi-theta);

DEAu = E1Au-E2Au;
DEAuErr = sqrt(E1AuErr^2+E2AuErr^2);
DEC = E2C-E1C;
DECErr = sqrt(E1CErr^2+E2CErr^2);

tC = 1/(a*d-b*c)*(d*DEAu-b*DEC)
tCerr = sqrt(1/(a*d-b*c)^2*(d^2*DEAuErr^2 +b^2*DECErr^2))
tAu = 1/(a*d-b*c)*(-c*DEAu+a*DEC)
tAuErr = 1/(a*d-b*c)*(c^2*DEAuErr^2+a^2*DECErr^2)

%% prøve 2

I = importfile('..\Data\Thickness\Sample2_0degree_480sec.asc');
E = 1:length(I);
Ierr = sqrt(I);
figure
peakBorders = [270,340;370,446];
name = 'samp2 pos 2'
data = fitGaussInSpectrum(E,I,Ierr,name,2,peakBorders);
E2C = c2E(data(1,1));
E2CErr = (data(2,1))/0.76535/1000;
E2Au = c2E(data(1,2));
E2AuErr = (data(2,2))/0.76535/1000;

I = importfile('..\Data\Thickness\Sample2_180degree_XXXsec.asc');
E = 1:length(I);
Ierr = sqrt(I);

peakBorders = [270,353;380,475];
name = 'samp2 pos 1'
data = fitGaussInSpectrum(E,I,Ierr,name,2,peakBorders);
legend('show')

E1C = c2E(data(1,1));
E1CErr = (data(2,1))/0.76535/1000;
E1Au = c2E(data(1,2));
E1AuErr = (data(2,2))/0.76535/1000;
%%

a = SC(E0)*K2(theta,mAu)+SC(E2Au)/cos(pi-theta);
b = 1/(2*cos(pi-theta))*(SAu(E2Au)-SAu(E1Au));
c = 1/(2*cos(pi-theta))*(SC(E1C)-SC(E2C));
d = SAu(E0)*K2(theta,mC)+SAu(E1C)/cos(pi-theta);

DEAu = E1Au-E2Au;
DEAuErr = sqrt(E1AuErr^2+E2AuErr^2);
DEC = E2C-E1C;
DECErr = sqrt(E1CErr^2+E2CErr^2);

tC = 1/(a*d-b*c)*(d*DEAu-b*DEC)
tCerr = sqrt(1/(a*d-b*c)^2*(d^2*DEAuErr^2 +b^2*DECErr^2))
tAu = 1/(a*d-b*c)*(-c*DEAu+a*DEC)
tAuErr = 1/(a*d-b*c)*(c^2*DEAuErr^2+a^2*DECErr^2)

function data = fitGaussInSpectrum(X,Y,Yerr,Name,n,peakBorder)
c2E = @(x) (x*0.76535+12.393);

% figure
hold on
xlim(c2E([200,500]))
ylim([0,max(Y(200:500))*1.1])

histogram('BinEdges',c2E([0,X(1:end)]),'BinCounts',Y,'EdgeColor','none','displayName',['Spec. ' Name])

xlabel('Energy (keV)')
ylabel('Counts (n)')
set(gca,'FontSize',14) 
% title(name)

name = {'C','Ag'};
shape = {'--','-.'};


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

beta0 = [(x2+x3)/2,(x3-x2)/3,max(y)/2,0,20];

w = 1./(yerr+1).^2;
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x',y,@fitfunction,beta0,'weights',w);

plot(c2E(x),fitfunction(beta,x),shape{i},'linewidth',2,'displayName',[name{i} ' fit'])
us = CovB/MSE;
mse =MSE;
MSECount(i) = MSE;
pValue(i) = 1-chi2cdf(MSE*(length(y)-5),(length(y)-5));
P_Value = 1-chi2cdf(MSE*(length(y)-5),(length(y)-5));

peakChannel(i) = beta(1,1);
peakUns(i) = us(1,1);
ci = nlparci(beta,R,'jacobian',J,'alpha',0.35);
peakUns(i) = (ci(1,2)-ci(1,1))/2;

end

peakUns;
peakChannel;
data = [peakChannel;peakUns;MSECount];
end


function y = fitfunction(beta,x)
    y = beta(4).*x+beta(5)+beta(3).*exp(-((x-beta(1))./(beta(2))).^2./2);
end
