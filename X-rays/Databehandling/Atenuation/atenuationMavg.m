clear all; close all;clc;
%%Kilder
% https://www.google.com/url?sa=i&source=images&cd=&ved=2ahUKEwi89LGLjPLlAhVCbFAKHTUkCmwQjRx6BAgBEAQ&url=https%3A%2F%2Fwww.researchgate.net%2Ffigure%2FA-Mass-attenuation-coefficients-of-various-heavy-metal-elements-and-x-ray-photon_fig1_327021329&psig=AOvVaw2a5nifCoMinO68vypk3eXY&ust=1574109250828331
%https://physics.nist.gov/PhysRefData/FFast/html/form.html
% Skalere I0 m�lingen med ca 6%
% Chi square fit 
%% Kalibrering 
a = 0.031179;
aUs = 3.605e-06;
b = -0.08267;
bUs = 0.0025276;
coVarAB = -8.2792e-09;
c2E = @(c) a*c+b;
c2EUs = @(c,cUs) sqrt((cUs*a).^2+(c*aUs).^2+bUs.^2+coVarAB.*c);
E2c = @(E) (E-b)/a


folderPath = '..\..\Absorbanter';
I0Name = '\Background\Background_36.5kV.txt';
I0Time = 73.74*(1-0.0586);

AbsorberNames ={{'\Absorber1\Absorber1  ca. 7 min.txt'},{'\Absorber2\Absorber2.txt'},...
    {'\Absorber3\Absorber3.txt','\Absorber3\Absorber3ca. 10 .txt'}};
AbsorberTimes = {[425.48],[599.75],[1532.94,600.14]};
AbsorberThickness = [33,30,25]*10^(-4); % Cm
AbsorberThicknessUs =[1,1,1]*10^(-4);



dat = load([folderPath I0Name]);
X = dat(:,1);
I=3:3:length(X);

E = c2E(X(I));
I0Y = dat(:,2);
I0Y = I0Y(I-2)+I0Y(I-1)+I0Y(I);

I0YUs = sqrt(I0Y)+(I0Y==0);
I0YT = I0Y/I0Time;
I0YTUs = I0YUs/I0Time;

%%
dat = load([folderPath '\GivenAtenuation\Gold.txt']);
AuE  = dat(:,1); % keV
AuAtt  = dat(:,2); %cm^-1
dat = load([folderPath '\GivenAtenuation\Molybdenum.txt']);
MoE  = dat(:,1); % keV
MoAtt  = dat(:,2); %cm^-1
dat = load([folderPath '\GivenAtenuation\Silver.txt']);
AgE  = dat(:,1); % keV
AgAtt  = dat(:,2); %cm^-1

%%
ESS = {MoE,AgE,AuE};
AttS = {MoAtt,AgAtt,AuAtt};

for i = 1:3
    ES = ESS{i};
    AttList = AttS{i};
    Att = interp1(ES,AttList,E);
    
    AbsorberName = AbsorberNames{i};
    AbsorberTime = sum(AbsorberTimes{i});
    AbsorberThicknes = AbsorberThickness(i);
    AbsorberThicknesUs = AbsorberThicknessUs(i);
    Y = zeros(size(X));
    for j =1:length(AbsorberName)
        dat = load([folderPath AbsorberName{j}]);
        Y = Y+dat(:,2);
    end
    Y = Y(I-2)+Y(I-1)+Y(I);

    YUs = sqrt(Y)+(Y==0);
    YT = Y/AbsorberTime;
    YTUs = YUs/AbsorberTime;   
    
    attenuation = (log(I0YT)-log(YT))/AbsorberThicknes;
    attenuationUs = sqrt(...
        (I0YTUs./(AbsorberThicknes*I0YT)).^2+ ...
        (YTUs./(AbsorberThicknes*YT)).^2+ ...
        ((log(I0YT)-log(YT))./AbsorberThicknes.^2*AbsorberThicknesUs).^2);
    
    %Plot af spectrum 
    figure
    hold on
    xlabel('Energy (E) [keV]')
    ylabel('Counts pr time [s^{-1}]')
    set(gca,'FontSize',15) 
    set(gca,'yscale','log')
    
    errorbar(E,I0YT,I0YTUs,'.','markersize',8)
    errorbar(E,YT,YTUs,'.','markersize',8)
   
    plot(E,I0YT.*exp(-Att.*AbsorberThicknes),'-b','linewidth',1)
    legend('I_0','I','Expected I','Location','northeast')

    %Plot attenuation
    figure
    xlabel('Energy (E) [keV]')
    ylabel('Attenuation coefficient \mu [\mu m^{-1}] ')
    set(gca,'FontSize',15) 
    hold on
    
    errorbar(E,attenuation,attenuationUs,'.','markersize',8)
    plot(ES,AttList,'-','linewidth',1)
    plot(E,Att,'-','linewidth',1)
    legend('Data','Expected','Location','northeast')

    %% Thikness
    th = (log(I0YT)-log(YT))./Att;
    thUs = sqrt(...
        (I0YTUs./(Att.*I0YT)).^2+ ...
        (YTUs./(Att.*YT)).^2);
    figure
    xlabel('Energy (E) [keV]')
    ylabel('Thikness [\mu m] ')
    set(gca,'FontSize',15) 
    hold on
    errorbar(E,th,thUs,'.','markersize',8)
    
    index = (E>8) & ~isnan(th) & (th~=Inf) & ~isnan(thUs)& (thUs~=Inf);
    
    W = 1./thUs(index).^2;
    
    % mu 
    T = sum(th(index).*W)./sum(W)*10000
    Tus = (sqrt(sum(W)))^-1*10000
    Tmean = mean(th(index))*10000
    Tstd = std(th(index))*10000
    figure
    histogram(th(index))
    
end