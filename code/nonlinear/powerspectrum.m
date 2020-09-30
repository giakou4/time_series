function [psM] = powerspectrum(xV,npoints,tittxt)
% [psM] = powerspectrum(xV,npoints,tittxt)
% POWERSPETRUM computes and plots the periodogram of a given
% a time series
% INPUTS:
%  xV      : vector of a scalar time series
%  npoints : number of samples to use (if larger than the length of xV, 
%            pad with zero)
%  tittxt  : string to be displayed in the title
% OUTPUT:
%  psM     : matrix of two columns: first column the frequencies, second
%            column the corresponding values of the periodogram.

n = length(xV);
if nargin == 1
    npoints = n;
    tittxt = '';
elseif nargin == 2
    tittxt = '';
end

[PxV,wV] = periodogram(xV,[],2*npoints);
fV = wV/(2*pi);
psV = 10*log10(PxV*pi); % One-sided, i.e. already multiplied by 2
figure(gcf)
clf
plot(fV,psV)
xlabel('frequency f')
ylabel('10*log_{10}(P_{per}(f))')
title([tittxt,' Periodogram (in dB), N=',int2str(npoints)])
ax = axis;
axis([0 0.5 ax(3) ax(4)])
psM = [fV';psV']';
