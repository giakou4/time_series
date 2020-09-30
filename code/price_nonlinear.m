%% GROUP
% Giakoumoglou Nikolaos AEM: 9043
% Gkrispanis Konstantinos AEM: 9306
%% Non-Linear analysis on Time-Series
%% Import data 
clc; clear all; close all;
data = importdata('ElectricPowerItaly.xls');
prices = data.data.prices(:,8);
iter = 1:24:length(prices);
prices = prices(iter); %Prices Northern Italy at 01.00 am
clear data; clear iter;
fig = 1;
Tmax = 1;
%% Previous work
% We also used a MA(2) filter to remove trend 
ma_filtre_trend = movingaveragesmooth(prices,2);
prices2 = prices - ma_filtre_trend;
prices2 = prices2(2:end-1);
% Best model was ARMA(1,2)
p = 1;
q = 2;
prices3 = prices2 - mean(prices2);
armamodel = armax(prices3,[p q]);
yfit = predict(armamodel,prices3,1);

figure(fig)
fig = fig + 1;
plot(prices3)
hold on
plot(yfit,'r')
legend('real data','fitted data')

residuals = prices3 - yfit;
n = length(residuals);

figure(fig)
fig = fig + 1;
plot(residuals)
title('residuals: real data substracted by fitted ARMA with 1-step prediction')
legend('residuals')
%% Create 20 iid time-series from <residuals>
IID = NaN*ones(n,21);
IID(:,1) = residuals;
for i = 2:21
    IID(:,i) = randsample(residuals,n,'false');
end
%% Linear features: Autocorrelation Function & Partial Autocorrelation Function 
taumax = 15;
ACF = NaN*ones(taumax+1,21);
PACF = NaN*ones(taumax,21);
for i = 1:21
    temp = autocorrelation(IID(:,i),taumax);
    ACF(:,i) = temp(:,2);
    
    PACF(:,i) = parautocor(IID(:,i),taumax);
end
close all;
%% Non-Linear features
%% Mutual Information
MI = NaN*ones(taumax+1,21);
for i = 1:21
    temp = mutualinformation(IID(:,i), taumax);
    MI(:,i) = temp(:,2);
end
close all;

figure(fig)
fig = fig + 1;
for i = 1:21
    subplot(5,5,i)
    plot(MI(:,i))
    title(sprintf('Mutual Information for lag %.0f',i));
    xlabel('\tau')
    ylabel('MI')
end

tau = 2; %for every plot we see lag                                                                                                                                                                                                                       = 2
%% False nearest neighbours
mmax = 10;
FNN = NaN*ones(mmax,21);
for i = 1:21
    temp = falsenearest(IID(:,i),tau,10,10,0);
    FNN(:,i) = temp(:,2);
end

figure(fig)
fig = fig + 1;
for i = 1:21
    subplot(5,5,i)
    plot(FNN(:,i))
    title(sprintf('Mutual Information for lag tau = %.0f',tau));
    xlabel('Embedding Dimension m')
    ylabel('False Nearest Neighbours')
end  

m = 2; %for every plot we see Embedding Dimension = 2
%% Correlation Dimension
V = NaN*ones(1,21);
for i = 1:21
      [~,~,~,~,nuM] = correlationdimension(IID(:,i),tau,mmax,'IID-1');
      V(1,i) = nuM(mmax,4);
end
%% Maximal Lyapunov Exponent
L = NaN*ones(1,21);
for i = 1:21
    temp = maxlyapunov(IID(:,i),tau,2,3,' ',20,5);
    L(1,i) = temp(1,1);
end
close all;
%% NRMSE for local fit
Tmax = 3;
NRMSE1 = NaN*ones(21,Tmax); %local linear: Tmax=1, nnei=1 q=0
for i = 1:21
    [NRMSE1(i,:),~] = localfitnrmse(IID(:,i),tau,m,Tmax,1,0);
end

NRMSE2 = NaN*ones(21,Tmax); %local linear: Tmax=1, nnei=2 q=0
for i = 1:21
    [NRMSE2(i,:),~] = localfitnrmse(IID(:,i),tau,m,Tmax,1,0);
end

NRMSE3 = NaN*ones(21,Tmax); %local linear: Tmax=1, nnei=2 q=3 -> OLS
for i = 1:21
    [NRMSE3(i,:),~] = localfitnrmse(IID(:,i),tau,m,Tmax,1,0);
end
%NRMSE all the same, keep 1st :)
NRMSE = NRMSE1;
%% Conclude:
% Linear: ACF, PACF
% Non-Linear: MI, FNN, V, L, NRMSE
%% ACF
figure(fig)
fig = fig+1;
suptitle('Histograms of Autocorrelation for lag=1..16')
for i=1:16
    figure(fig-1)
    subplot(4,4,i)
    temp = ACF(i,:);
    histogram(temp)
    title(sprintf('lag=%.0f',i))
    hold on
    plot([ACF(i,1) ACF(i,1)],[0 10],'r','LineWidth',4)
    legend('20 IID from residuals','residuals')
end
%% PACF
figure(fig)
fig = fig+1;
suptitle('Histograms of Partial Autocorrelation for lag=1..15')
for i=1:15
    figure(fig-1)
    subplot(4,4,i)
    temp = PACF(i,:);
    histogram(temp)
    title(sprintf('lag=%.0f',i))
    hold on
    plot([PACF(i,1) PACF(i,1)],[0 10],'r','LineWidth',4)
    legend('20 IID from residuals','residuals')
end
%% MI
figure(fig)
fig = fig+1;
suptitle('Histograms of Mutual Information for lag=1..16')
for i=1:16
    figure(fig-1)
    subplot(4,4,i)
    temp = MI(i,:);
    histogram(temp)
    title(sprintf('lag=%.0f',i))
    hold on
    plot([MI(i,1) MI(i,1)],[0 10],'r','LineWidth',4)
    legend('20 IID from residuals','residuals')
end
%% FNN
figure(fig)
fig = fig+1;
suptitle('Histograms of False Nearest Neighbours for lag=1..16')
for i=1:3
    figure(fig-1)
    subplot(1,3,i)
    temp = FNN(i,:);
    histogram(temp)
    title(sprintf('lag=%.0f',i))
    hold on
    plot([FNN(i,1) FNN(i,1)],[0 10],'r','LineWidth',4)
    legend('20 IID from residuals','residuals')
end
%% V
figure(fig)
fig = fig+1;
histogram(V(2:end))
hold on
plot([V(1) V(1)],[0 12],'r','LineWidth',4)
legend('20 IID from residuals','residuals')
title('Histogram of estimated Correlation Dimension')
%% L
figure(fig)
fig = fig+1;
histogram(L(2:end))
hold on
plot([L(1) L(1)],[0 12],'r','LineWidth',4)
legend('20 IID from residuals','residuals')
title('Histogram of estimated Maximal Lyapunov Exponents')
%% NRMSE
figure(fig)
fig = fig+1;
suptitle('Histograms of NRMSE for Tmax=1...3')
for i=1:3
    figure(fig-1)
    subplot(3,1,i)
    hist(NRMSE(2:end,i))
    hold on 
    plot([NRMSE(1,i) NRMSE(1,i)],[0 12],'r','LineWidth',4)
    legend('20 IID from residuals','residuals')
    title(sprintf('Tmax = %.0f',i))
end