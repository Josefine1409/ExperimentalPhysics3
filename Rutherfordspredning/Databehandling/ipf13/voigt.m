function g=voigt(x,pos,Gausswidth,LorentzWidth)
% Unit height Voigt profile function. x is the independent variable
% (energy, wavelength, etc), Gausswidth is the Doppler (Gaussian) width,
% and LorentzWidth is the pressure (LorentzWidth) width. 
% Version 2, June 13, 2018
%
% Example:
% x=0:.01:2;
% pos=1.2;
% Gausswidth=.2;
% LorentzWidth=.1;
% v=voigt(x,pos,Gausswidth,LorentzWidth);
% plot(x,gaussian(x,pos,Gausswidth),x,lorentzian(x,pos,LorentzWidth),x,v)
%
gau=gaussian(x-pos,0,Gausswidth);
lor=lorentzian(x,pos,LorentzWidth);
VoigtConv=conv(gau,lor,'full');
g=VoigtConv./max(VoigtConv);
g=g(1:2:length(g));

function g = gaussian(x,pos,wid)
%  gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 1988
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.60056120439323.*wid)) .^2);

function g = lorentzian(x,position,width)
% lorentzian(x,position,width) Lorentzian function.
% where x may be scalar, vector, or matrix
% position and width scalar
% T. C. O'Haver, 1988
% Example: lorentzian([1 2 3],2,2) gives result [0.5 1 0.5]
g=ones(size(x))./(1+((x-position)./(0.5.*width)).^2);
