% The program determines an AR process, by determining the roots of the
% characteristic polynomial, and computes (and plots if specified) for the
% given AR coefficients the autocorrelation and partial autocorrelation.
clear all
display = 1; % if 1 then display graph
maxtau = 20; % maximum lag to show 
% The roots of the characteristic polynomial of AR(p)
% B^p + phi_1 * B^(p-1) + ... + phi_(p-1) * B + phi_p = 0
% rooV = [0.5+0.8*i 0.5-0.8*i -0.5+0.8*i -0.5-0.8*i]'; 
rooV = [0.5+0.8*i 0.5-0.8*i 1]'; 
p = length(rooV); 
tmpV = poly(rooV);
phiV = -tmpV(2:p+1)'; % The phi coefficients of AR(p)

rhoV = ARcoef2autocorr(phiV,maxtau,display); % The autocorrelation from phi
pacfV = acf2pacf(rhoV,display); % The partial autocorrelation from autocorrelation


