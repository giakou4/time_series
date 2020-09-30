function pautV = parautocor(xV,maxtau)
% pautV = parautocor(xV,maxtau)
% PARAUTOCOR computes the partial autocorrelation on the time
% series 'xV' for delays up to a maximum lag 'maxtau'. Autoregressive
% models of orders 1,2,...,maxtau are fitted to the time series and the
% last coefficient for each model is the respective partial
% autocorrelation, for lags 1,2,...,maxtau.
% INPUT
% - xV      : a vector for the time series
% - maxtau  : the maximum lag to compute partial autocorrelation for
% OUTPUT
% - pauV    : the vector of the partial autocorrelations for the given lags

n = length(xV);
if n<maxtau+1
    return;
end
xM = NaN*ones(n-maxtau,maxtau);
for i=1:maxtau
    xM(i+1:n,i)=xV(1:n-i);
end
pautV = NaN*ones(maxtau,1);
for tau=1:maxtau
    [Q,R]=qr([ones(n-tau,1) xM(tau+1:n,1:tau)],0);
    phiV =  R\(Q'*xV(tau+1:n));
    pautV(tau) = phiV(end);
end
