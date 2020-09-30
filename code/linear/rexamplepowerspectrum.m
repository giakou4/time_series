% Simulation of a SARMA process. 
% 1. Generation of a time series from a given SARMA process
% 2. Autocorrelation and partial autocorrelation
% 3. Estimation of power spectrum: 
%    a) classical approach (periodogram), 
%    b) parametric approach (AR model). clear all
rooV = [0.4+0.3*i 0.4-0.3*i]';
phi0= 0;
thetaV = [-5 0.4];
phisV = [0.8];
thetasV = [];
s = 6;
n = 1000;
sdnoise = 1;
maxtau = 100;
Tmax = 10;
proptest = 0.3;
alpha = 0.05;

tmpV = poly(rooV);
phiV = [phi0; -tmpV(2:length(rooV)+1)'];
pgen = length(phiV)-1;
qgen = length(thetaV);
psgen = length(phisV);
qsgen = length(thetasV); 
% 1. Generation of a time series from a given SARMA process
[xV,phiallV,thetaallV]= generateSARMAts(phiV,thetaV,phisV,thetasV,s,n,sdnoise);
figno = 1;
figure(figno)
clf
plot(xV,'.-')
hold on
xlabel('t')
ylabel('x(t)')
title(sprintf('SARMA(%d,%d)x(%d,%d)_%d, time history',pgen,qgen,psgen,qsgen,s))

% 2. Autocorrelation and partial autocorrelation 
% 2a. Autocorrelation
[acM] = autocorrelation(xV, maxtau);
zalpha = norminv(1-alpha/2);
autlim = zalpha/sqrt(n);
figno = figno + 1;
figure(figno)
clf
hold on
for ii=1:maxtau
    plot(acM(ii+1,1)*[1 1],[0 acM(ii+1,2)],'b','linewidth',1.5)
end
plot([0 maxtau+1],[0 0],'k','linewidth',1.5)
plot([0 maxtau+1],autlim*[1 1],'--c','linewidth',1.5)
plot([0 maxtau+1],-autlim*[1 1],'--c','linewidth',1.5)
xlabel('\tau')
ylabel('r(\tau)')
title(sprintf('SARMA(%d,%d)x(%d,%d)_%d, autocorrelation',pgen,qgen,psgen,qsgen,s))

% 2b. Partial autocorrelation
display = 1;
pacfV = parautocor(xV,maxtau);
figno = figno + 1;
figure(figno)
clf
hold on
for ii=1:maxtau
    plot(acM(ii+1,1)*[1 1],[0 pacfV(ii)],'b','linewidth',1.5)
end
plot([0 maxtau+1],[0 0],'k','linewidth',1.5)
plot([0 maxtau+1],autlim*[1 1],'--c','linewidth',1.5)
plot([0 maxtau+1],-autlim*[1 1],'--c','linewidth',1.5)
xlabel('\tau')
ylabel('\phi_{\tau,\tau}')
title(sprintf('SARMA(%d,%d)x(%d,%d)_%d, partial autocorrelation',pgen,qgen,psgen,qsgen,s))

% 3. Estimation of power spectrum: 
%    a) classical approach (periodogram), 
[PxV,wV] = periodogram(xV);
fV = wV/(2*pi);
PperV = 10*log10(PxV*pi); % One-sided, 
figno = figno + 1;
figure(figno)
clf
plot(fV,PperV)
hold on
xlabel('frequency f')
ylabel('10*log10(P(f))')
title(sprintf('SARMA(%d,%d)x(%d,%d)_%d, power spectral density estimation',pgen,qgen,psgen,qsgen,s))
% axis([0 0.5 -30 35])
% ax = axis;

%    b) parametric approach (AR model)
p = input('Give the order p of the AR model for power spectral density estimation >'); 
[P2xV,w2V] = pmcov(xV,p);
f2V = w2V/(2*pi);
ParV = 10*log10(P2xV*pi); % One-sided, 
figure(figno)
plot(f2V,ParV,'r')
legend('P_{per}',sprintf('P_{AR(%d)}',p),'Location','Best');

