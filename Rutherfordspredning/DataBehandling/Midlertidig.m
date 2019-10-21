Ein = 349.9

ts=[300 300 300 300 350 600 600 600 600 600 600 600];
FCs=[46449 51781 60300 65892 81962 34228 103355 80585 101501 106019 102580 37053];

GCs=CountsAg;
CCs=CountsC;
GEs=EnergyAg;
CEs=EnergyC;

sigmaGCs=CountsUnsAg;
sigmaCCs=CountsUnsC;
sigmaGEs=EnergyUnsAg;
sigmaCEs=EnergyUnsC;

GCs=GCs./FCs;
CCs=CCs./FCs;
sigmaGCs = sigmaGCs./FCs;
sigmaCCs = sigmaCCs./FCs;

GCs=GCs([1:5,7:11])
CCs=CCs([1:5,7:11])
sigmaGCs = sigmaGCs([1:5,7:11])
sigmaCCs = sigmaCCs([1:5,7:11])


cs=@(theta,m,E) 1./(sin(theta/2)).^4.*1./(K2(theta,m).*E).^2;

% [betaG,RG,JG,CovB,MSE]=nlinfit(thetas,GCs,@(C,thetas) C(1)*cs(thetas,mG,Ein)+C(2),[1,0]);
% [betaC,RC,JC,CovB,MSE]=nlinfit(thetas,CCs,@(C,thetas) C(1)*cs(thetas,mC,Ein)+C(2),[1,0]);
[betaG,RG,JG,CovB,MSE]=nlinfit(thetas([1:5,7:11]),GCs,@(C,thetas) C(1)*cs(thetas,mG,Ein)+C(2),[1,0]);
[betaC,RC,JC,CovB,MSE]=nlinfit(thetas([1:5,7:11]),CCs,@(C,thetas) C(1)*cs(thetas,mC,Ein)+C(2),[1,0]);

Thetas=linspace(thetas(1),thetas(end),1000);

figure
errorbar(theta,GEs,sigmaGEs,'r.','markersize',10)
hold on
xlabel('Scattering Angle [deg]')
ylabel('Scattering Energy [KeV]')
title('Scattering on Gold')
plot(Thetas/pi*180,EoutNy(Ein,Thetas,1),'linewidth',2)

figure
errorbar(theta,CEs,sigmaCEs,'r.','markersize',10)
hold on
xlabel('Scattering Angle [deg]')
ylabel('Scattering Energy [KeV]')
title('Scattering on Carbon')
plot(Thetas/pi*180,EoutNy(Ein,Thetas,2),'linewidth',2)

% 
% figure
% [Ypred,deltaY] = nlpredci(@(C,th) C(1)*cs(th,mG,Ein)+C(2),Thetas,betaG,RG,'jacobian',JG,'alpha',0.35);
% 
% errorbar(thetas,GCs,sigmaGCs,'r.','markersize',10)
% hold on
% xlabel('Scattering Angle')
% ylabel('Normed Counts')
% title('Scattering on Gold')
% plot(Thetas,Ypred,'linewidth',2)
% 
% plot(Thetas,Ypred+deltaY,'k--','linewidth',1)
% plot(Thetas,Ypred-deltaY,'k--','linewidth',1)
% legend('Fit','Data','Fit confidence','Location','southwest')
% 
% 
% figure
% [Ypred,deltaY] = nlpredci(@(C,th) C(1)*cs(th,mC,Ein)+C(2),Thetas,betaC,RC,'jacobian',JC,'alpha',0.35);
% 
% errorbar(thetas,CCs,sigmaCCs,'r.','markersize',10)
% hold on
% xlabel('Scattering Angle [deg]')
% ylabel('Normed Counts')
% title('Scattering on Carbon')
% plot(Thetas,Ypred,'linewidth',2)
% 
% plot(Thetas,Ypred+deltaY,'k--','linewidth',1)
% plot(Thetas,Ypred-deltaY,'k--','linewidth',1)
% legend('Fit','Data','Fit confidence','Location','southwest')

figure
[Ypred,deltaY] = nlpredci(@(C,th) C(1)*cs(th,mG,Ein)+C(2),Thetas,betaG,RG,'jacobian',JG,'alpha',0.35);

hold on
xlabel('Scattering Angle')
ylabel('Normed Counts')
title('Scattering on Gold')
plot(Thetas/pi*180,Ypred,'linewidth',2)
errorbar(theta([1:5,7:11]),GCs,sigmaGCs,'r.','markersize',10)

plot(Thetas/pi*180,Ypred+deltaY,'k--','linewidth',1)
plot(Thetas/pi*180,Ypred-deltaY,'k--','linewidth',1)
legend('Fit','Data','Fit confidence','Location','southwest')


figure
[Ypred,deltaY] = nlpredci(@(C,th) C(1)*cs(th,mC,Ein)+C(2),Thetas,betaC,RC,'jacobian',JC,'alpha',0.35);

hold on
xlabel('Scattering Angle [deg]')
ylabel('Normed Counts')
title('Scattering on Carbon')
plot(Thetas/pi*180,Ypred,'linewidth',2)
errorbar(theta([1:5,7:11]),CCs,sigmaCCs,'r.','markersize',10)

plot(Thetas/pi*180,Ypred+deltaY,'k--','linewidth',1)
plot(Thetas/pi*180,Ypred-deltaY,'k--','linewidth',1)


legend('Fit','Data','Fit confidence','Location','southwest')


%% Fit
function eout = EoutNy(Ein,Theta,id)
Ein = Ein/1000
mG=196.966-4.4858e-4*79-0.03343120468;
mC=12-4.4858e-4*12;
m1=1.007276;

TAu = 25e-10;
TC = 250e-10;
eout =[];
m = [mG,mC];
K2=@(theta)((m1*cos(theta)+sqrt(m(id)^2-m1^2*sin(theta).^2))./(m1+m(id))).^2;
for theta = Theta;
SC = stoppingpowerC(Ein);
SAu = stoppingpowerAu(Ein);
if (id==1)
    Escatter = (Ein-TAu/2*SAu).*K2(theta);
    SC = stoppingpowerC(Escatter);
    SAu = stoppingpowerAu(Escatter);
    if (theta<pi/2)
        eout(end+1) = Escatter-(TAu./2.*SAu+TC.*SC)./norm(cos(theta));
    else
        eout(end+1) = Escatter-(TAu./2.*SAu)./norm(cos(theta));
    end
else
    Escatter = (Ein-TAu*SAu-TC*SC).*K2(theta);
    SC = stoppingpowerC(Escatter);
    SAu = stoppingpowerAu(Escatter);
    if (theta<pi/2)
        eout(end+1) = Escatter-(TC/2.*SC)./norm(cos(theta));
    else
        eout(end+1) = Escatter-(TC/2.*SC+TAu.*SAu)./norm(cos(theta));
    end
end
end
eout = eout*1000;
end