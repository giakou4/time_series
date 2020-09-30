% Simulation of an ARMA process. 
% 1. Generation of a time series from a given ARMA process
% 2. Autocorrelation and partial autocorrelation 
% 3. Fit of a ARMA model.
% 4. Prediction error statistic for the ARMA model.
% 5. Multi-step predictions from a given target point with the ARMA model.
clear all
rooV = [0.9]';
phi0= 0;
thetaV = [];
n = 100;
sdnoise = 1;
maxtau = 10;
Tmax = 10;
proptest = 0.3;
alpha = 0.05;

tmpV = poly(rooV);
phiV = [phi0; -tmpV(2:length(rooV)+1)'];
pgen = length(phiV)-1;
qgen = length(thetaV);
% 1. Generation of a time series from a given ARMA process
xV = generateARMAts(phiV,thetaV,n,sdnoise);
figno = 1;
figure(figno)
clf
plot(xV,'.-')
hold on
xlabel('t')
ylabel('x(t)')
title(sprintf('ARMA(%d,%d), time history',pgen,qgen))

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
title(sprintf('ARMA(%d,%d), autocorrelation',pgen,qgen))

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
title(sprintf('ARMA(%d,%d), partial autocorrelation',pgen,qgen))

% 3. Fit of an ARMA model.
p = input('Give the order p of the AR part >'); 
q = input('Give the order q of the MA part >'); 
[nrmseV,phiallV,thetaallV,SDz,aicS,fpeS]=fitARMA(xV,p,q,Tmax);
fprintf('===== ARMA model ===== \n');
fprintf('Estimated coefficients of phi(B):\n');
disp(phiallV')
fprintf('Estimated coefficients of theta(B):\n');
disp(thetaallV')
fprintf('SD of noise: %f \n',SDz);
fprintf('AIC: %f \n',aicS);
fprintf('FPE: %f \n',fpeS);
fprintf('\t T \t\t NRMSE \n');
disp([[1:Tmax]' nrmseV])
if Tmax>3
    figno = figno + 1;
    figure(figno)
    clf
    plot([1:Tmax]',nrmseV,'.-k')
    hold on
    plot([1 Tmax],[1 1],'y')
    xlabel('T')
    ylabel('NRMSE')
    title(sprintf('ARMA(%d,%d), fitting error',p,q))
end    

% 4. Prediction error statistic for the ARMA model.
nlast = proptest*n;
tittxt = sprintf('ARMA(%d,%d), %%test=%1.2f, prediction error',p,q,proptest);
figno = figno + 1;
figure(figno);
clf
[nrmseV,preM] = predictARMAnrmse(xV,p,q,Tmax,nlast,'example');
figno = figno + 1;
figure(figno);
clf
plot([n-nlast+1:n]',xV(n-nlast+1:n),'.-')
hold on
plot([n-nlast+1:n]',preM(:,1),'.-r')
if Tmax>1
    plot([n-nlast+1:n]',preM(:,2),'.-c')
	if Tmax>2
        plot([n-nlast+1:n]',preM(:,3),'.-k')
    end
end
switch Tmax
    case 1
        legend('true','T=1','Location','Best')
    case 2
        legend('true','T=1','T=2','Location','Best')
    otherwise
        legend('true','T=1','T=2','T=3','Location','Best')
end
% 5. Multi-step predictions from a given target point with the ARMA model.
n1 = n-Tmax;
figno = figno + 1;
figure(figno);
clf
[preV] = predictARMAmultistep(xV,n1,p,q,Tmax,'example');
