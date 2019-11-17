clear all; close all;clc;
%% Kalibrering 
a = 0.031179;
aUs = 3.605e-06;
b = -0.08267;
bUs = 0.0025276;

c2E = @(c) a*c+b;
c2EUs = @(c,cUs) sqrt((cUs*a).^2+(c*aUs).^2+bUs.^2);
E2c = @(E) (E-b)/a


folderPath = '..\..\Absorbanter';
I0Name = '\Background\Background_36.5kV.txt';
I0Time = 73.74;

AbsorberNames ={{'\Absorber1\Absorber1  ca. 7 min.txt'},{'\Absorber2\Absorber2.txt'},...
    {'\Absorber3\Absorber3.txt','\Absorber3\Absorber3ca. 10 .txt'}};
AbsorberTimes = {[425.48],[599.75],[1532.94,600.14]};
AbsorberThickness = [1,1,1];
AbsorberThicknessUs =[0,0,0];


dat = load([folderPath I0Name]);
X = dat(:,1);
E = c2E(X);
I0Y = dat(:,2);
I0YUs = sqrt(I0Y)+(I0Y==0);
I0YT = I0Y/I0Time;
I0YTUs = I0YUs/I0Time;


for i = 1:3
    AbsorberName = AbsorberNames{i};
    AbsorberTime = sum(AbsorberTimes{i});
    AbsorberThicknes = AbsorberThickness(i);
    AbsorberThicknesUs = AbsorberThicknessUs(i);
    Y = zeros(size(E));
    for j =1:length(AbsorberName)
        dat = load([folderPath AbsorberName{j}]);
        Y = Y+dat(:,2);
    end
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
    legend('I_0','I','Location','northeast')
    
    %Plot e
    figure
    xlabel('Energy (E) [keV]')
    ylabel('attenuation coefficient [\mu m^{-1}] ')
    set(gca,'FontSize',15) 
    hold on
    errorbar(E,attenuation,attenuationUs,'.','markersize',8)
end