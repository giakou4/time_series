% The program determines an MA process, by determining the roots of the
% characteristic polynomial, and computes (and plots if specified) for the
% given MA coefficients the autocorrelation and partial autocorrelation.
clear all
display = 1; % if 1 then display graph
maxtau = 20; % maximum lag to show 
% The roots of the characteristic polynomial of MA(q)
% B^p - theta_1 * B^(q-1) + ... + theta_(q-1) * B + theta_q = 0
rooV = [0.5+0.5*i 0.5-0.5*i -0.3+0.5*i -0.3-0.5*i]';
q = length(rooV);
tmpV = poly(rooV);
thetaV = -tmpV(2:length(rooV)+1)'; % The theta coefficients of MA(q)

rhoV = zeros(maxtau,1);
rhoV(1:q) = MAcoef2autocorr(thetaV,display); % The autocorrelation from theta
pacfV = acf2pacf(rhoV,display); % The partial autocorrelation from autocorrelation


