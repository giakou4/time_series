%% GROUP
% Giakoumoglou Nikolaos AEM: 9043
% Gkrispanis Konstantinos AEM: 9306
%% Linear analysis on Time-Series: Demand
%% Import data 
clc; close all; clear all;
data = importdata('ElectricPowerItaly.xls');
demand = data.data.demand(:,8);
iter = 1:24:length(demand);
demand = demand(iter); %Demand Northern Italy at 01.00 am
clear data; clear iter;
n = length(demand);
fig = 1;
alpha = 0.05;
maxtau = 100; % for ACF and PACF
Tmax = 1; % predict 1 day ahead
%% Demand: raw data
figure(fig)
fig = fig+1;
plot(demand)
title('Demand Northern Italy at 01.00am raw data [9043-9306]')
xlabel('hour')
ylabel('demand')
[H,P] = adftest(demand);
disp(['Demand: raw data are stationary if H=1 ---> H=',num2str(H)])
%% Autocorrelation Function (ACF)
[ACF] = autocorrelation(demand, maxtau);
zalpha = norminv(1-alpha/2);
autlim = zalpha/sqrt(n);
figure(fig)
fig = fig + 1;
clf
hold on
for ii=1:maxtau
    plot(ACF(ii+1,1)*[1 1],[0 ACF(ii+1,2)],'b','linewidth',1.5)
end
plot([0 maxtau+1],[0 0],'k','linewidth',1.5)
plot([0 maxtau+1],autlim*[1 1],'--c','linewidth',1.5)
plot([0 maxtau+1],-autlim*[1 1],'--c','linewidth',1.5)
xlabel('\tau')
ylabel('r(\tau)')
title('Autocorrelation of demand: raw data[9043-9306]')
%% Partial Autocorrelation Function (PACF)
PACF = parautocor(demand,maxtau);
figure(fig)
fig = fig + 1;
clf
hold on
for ii=1:maxtau
    plot(ACF(ii+1,1)*[1 1],[0 PACF(ii)],'b','linewidth',1.5)
end
plot([0 maxtau+1],[0 0],'k','linewidth',1.5)
plot([0 maxtau+1],autlim*[1 1],'--c','linewidth',1.5)
plot([0 maxtau+1],-autlim*[1 1],'--c','linewidth',1.5)
xlabel('\tau')
ylabel('\phi_{\tau,\tau}')
title('Partial Autocorrelation of demand: raw data [9043-9306]')
%% Remove trend
ma_filtre_trend = movingaveragesmooth(demand,2);
demand2 = demand - ma_filtre_trend;
figure(fig)
fig = fig+1;
plot(demand2);
demand2 = demand2(2:end-1);
title('Demand, trend is removed [9043-9306]')
xlabel('hours')
ylabel('demand')
[H,P] = adftest(demand2);
disp(['Demand after detrending are stationary if H=1 ---> H=',num2str(H)])
%% Autocorrelation Function (ACF)
[ACF] = autocorrelation(demand2, maxtau);
zalpha = norminv(1-alpha/2);
autlim = zalpha/sqrt(n);
figure(fig)
fig = fig + 1;
clf
hold on
for ii=1:maxtau
    plot(ACF(ii+1,1)*[1 1],[0 ACF(ii+1,2)],'b','linewidth',1.5)
end
plot([0 maxtau+1],[0 0],'k','linewidth',1.5)
plot([0 maxtau+1],autlim*[1 1],'--c','linewidth',1.5)
plot([0 maxtau+1],-autlim*[1 1],'--c','linewidth',1.5)
xlabel('\tau')
ylabel('r(\tau)')
title('Autocorrelation of demand, trend is removed [9043-9306]')
%% Partial Autocorrelation Function (PACF)
PACF = parautocor(demand2,maxtau);
figure(fig)
fig = fig + 1;
clf
hold on
for ii=1:maxtau
    plot(ACF(ii+1,1)*[1 1],[0 PACF(ii)],'b','linewidth',1.5)
end
plot([0 maxtau+1],[0 0],'k','linewidth',1.5)
plot([0 maxtau+1],autlim*[1 1],'--c','linewidth',1.5)
plot([0 maxtau+1],-autlim*[1 1],'--c','linewidth',1.5)
xlabel('\tau')
ylabel('\phi_{\tau,\tau}')
title('Partial Autocorrelation of demand, trend is removed [9043-9306]')
%% Remove seasonality
ma_filtre_seasonality = movingaverageseasonal(demand2,7);
demand3 = demand2 - ma_filtre_seasonality;

figure(fig)
fig = fig+1;
plot(demand3);
title('Demand, seasonality is removed [9043-9306]')
xlabel('hours')
ylabel('demand')
demand3 = demand3(4:end-3);
[H,P] = adftest(demand3);
disp(['Demand after detrending are removing seasonality are stationary if H=1 ---> H=',num2str(H)])
%% Autocorrelation Function (ACF)
[ACF] = autocorrelation(demand3, maxtau);
zalpha = norminv(1-alpha/2);
autlim = zalpha/sqrt(n);
figure(fig)
fig = fig + 1;
clf
hold on
for ii=1:maxtau
    plot(ACF(ii+1,1)*[1 1],[0 ACF(ii+1,2)],'b','linewidth',1.5)
end
plot([0 maxtau+1],[0 0],'k','linewidth',1.5)
plot([0 maxtau+1],autlim*[1 1],'--c','linewidth',1.5)
plot([0 maxtau+1],-autlim*[1 1],'--c','linewidth',1.5)
xlabel('\tau')
ylabel('r(\tau)')
title('Autocorrelation of demand, trend and seasonality are removed [9043-9306]')
%% Partial Autocorrelation Function (PACF)
PACF = parautocor(demand3,maxtau);
figure(fig)
fig = fig + 1;
clf
hold on
for ii=1:maxtau
    plot(ACF(ii+1,1)*[1 1],[0 PACF(ii)],'b','linewidth',1.5)
end
plot([0 maxtau+1],[0 0],'k','linewidth',1.5)
plot([0 maxtau+1],autlim*[1 1],'--c','linewidth',1.5)
plot([0 maxtau+1],-autlim*[1 1],'--c','linewidth',1.5)
xlabel('\tau')
ylabel('\phi_{\tau,\tau}')
title('Partial Autocorrelation of demand, trend and seasonality are removed [9043-9306]')
%% TILL HERE:
% Trend and Seasonality are alreaedy been removed from the time series
% <demand3> is NOT white noise at a=0.05 significance
% <demand3> holds the stationary time-series
% <ma_filtre_trend> holds the trend
% <ma_filtre_seasonality> holds the seasonality
%% Fit ARMA model using Brute-Force by finding minimum AIC
brute_force = NaN*ones(36,3);
i = 1;
for p=0:10
    for q=0:10
        [nrmseV,phiallV,thetaallV,SDz,aicS,fpeS]=fitARMA(demand3,p,q,Tmax);
        brute_force(i,1) = p;
        brute_force(i,2) = q;
        brute_force(i,3) = aicS;
        i = i+1;
    end
end

[minAIC index] = min(brute_force(:,3));
%p = brute_force(index,1);
%q = brute_force(index,2);
fprintf('Minimum AIC after Brute-Force appears for p=%.0f q=%.0f\n',p,q);

figure(fig)
fig = fig+1
for ii=1:10:(121-10)
    hold on
    plot(1:10,brute_force(ii:ii+9,3))
end
ylabel('AIC')
xlabel('q')
legend('p=1','p=2','p=3','p=4','p=5','p=6','p=7','p=8','p=9','p=10')
title('AIC for all combination p=0:10 q=0:10 [9043-9306]')
%% Fit my ARMA model
p = 4;
q = 3;
[nrmseV,phiallV,thetaallV,SDz,aicS,fpeS]=fitARMA(demand3,p,q,Tmax);
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
%% Find a more economical ARMA model
%{
stdBF = std(brute_force(:,3));
sumpq = p+q; 
econ_p = -1;
econ_q = -1;
for i = 1:121
   if( brute_force(i,3)<(minAIC + stdBF/8) & brute_force(i,3)>(minAIC - stdBF/8) )
       if(brute_force(i,1)+brute_force(i,2)<sumpq)
           econ_p = brute_force(i,1);
           econ_q = brute_force(i,2);
           sumpq = econ_p + econ_q;
       end
   end   
end
fprintf('Economy ARMA: p = %.0f q = %.0f\n',econ_p,econ_q);
clear stdBF; clear sumpq; clear i;
%}
%% 1 step prediction with my ARMA model
proptest = 0.90;
n = length(demand3);
nlast = floor(proptest*n);
[nrmseV,preM] = predictARMAnrmse(demand3,p,q,Tmax,nlast);
figure(fig);
fig = fig + 1;
clf
plot([n-nlast+1:n]',demand3(n-nlast+1:n),'.-')
hold on
plot([n-nlast+1:n]',preM(:,1),'.-r')
legend('true','T=1','Location','Best')
title(sprintf('Fit ARMA(%.0f,%.0f) on stationary demand time-series and predict 1 step ahead [9043-9306]',p,q));
%% Fit AR(5)
p = 5;
q = 0;
[nrmseV,phiallV,thetaallV,SDz,aicS,fpeS]=fitARMA(demand3,p,q,Tmax);
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
%% Fit AR(5) model & 1 step prediction with AR(5) for proptest = 0.05 to 0.8 with step=0.05
n = length(demand3);
proptest = 0.05:0.05:0.8;
NRMSE_array = ones(length(proptest),1);
for i = 1:length(proptest);
    nlast = floor(proptest(i)*n);    
    [nrmseV,preM] = predictARMAnrmse(demand3,p,q,Tmax,nlast);
    NRMSE_array(i) = nrmseV;
    figure(fig);
    subplot(4,4,i);
    plot([n-nlast+1:n]',demand3(n-nlast+1:n),'.-')
    hold on
    plot([n-nlast+1:n]',preM(:,1),'.-r')
    legend('true','T=1','Location','Best')
    title(sprintf('Fit ARMA(%.0f,%.0f) on stationary demand time-series and predict 1 step ahead\ntest-set is the last %.2f%% of observations\nNRMSE=%f [9043-9306]',p,q,100*proptest(i),nrmseV));
end
fig = fig + 1;
figure(fig)
plot(100.*proptest,NRMSE_array)
title('NRMSE vs test set')
xlabel('size of test set: % of total size time-series')
ylabel('NRMSE')
