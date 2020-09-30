function [nrmseV,preM,phiV] = linearpredictnrmse(xV,m,Tmax,nlast,tittxt)
% [nrmseV,phiV] = linearpredictnrmse(xV,m,Tmax,nlast,tittxt)
% LINEARPREDICTNRMSE make predictions with an AR model on a last part
% of a given time series and computes the prediction error (NRMSE measure)
% for T-step ahead predictions.
% INPUTS:
%  xV      : vector of the scalar time series
%  m       : the embedding dimension.
%  Tmax    : the prediction horizon, the predict is made for T=1...Tmax steps
%            ahead.
%  nlast   : the size of the test set to compute the prediction error on
%          : if not specified, it is half the length of the time series
%  tittxt  : string to be displayed in the title of the figure 
%            if not specified, no plot is made
% OUTPUT: 
%  nrmseV  : vector of length Tmax, the nrmse for the predictions for T-mappings, 
%            T=1...Tmax, on the test set
%  preM    : matrix of nlast columns and Tmax rows, having the T-ahead 
%            predictions at column T, T=1...Tmax (the first T-1 components
%            are NaN).
%  phiV    : the coefficients of the estimated AR time series (of length (m+1)
%            with phi(0) as first component

sizeofmark = 15; 
n = length(xV);
if nargin==4
    tittxt = [];
elseif nargin==3
    tittxt = [];
    nlast = round(n/2);
end
if nlast>=n-2*m,
    error('test set is too large for the given time series!')
end
n1 = n-nlast;  % size of training set
x1V = xV(1:n1); 
mx1 = mean(x1V(1:n1)); 
y1V = x1V(1:n1)-mx1; % centralize to 0 the training set
% th = ar(y1V,m);  % This is used in a previous version!
% aV = get(th,'a');
aV = armcov(y1V,m);
aV = -aV(2:m+1);
a0 = (1-sum(aV))*mx1; 
phiV = [a0 aV]';

preM = NaN*ones(n+Tmax-1,Tmax); % for simplicity use the indices for the whole
                                % time series, the first n1 will be ignored
for i=n1:n-1
    preV = NaN*ones(m+Tmax,1);
    preV(1:m)=xV(i-m+1:i)-mx1;
    for T=1:Tmax
        preV(m+T)=aV*preV(m+T-1:-1:T);
        preM(i+T,T)=preV(m+T);    
    end
end
preM = preM + mx1*ones(size(preM));
nrmseV = ones(Tmax,1);
for T=1:Tmax
    nrmseV(T) = nrmse(xV(n1+T:n),preM(n1+T:n,T));
end
preM = preM(n1+1:n,:);

if ~isempty(tittxt)
	figno = gcf;
	figure(figno)
	clf
	plot([1:Tmax]',nrmseV,'k')
	hold on
	plot([1:Tmax]',nrmseV,'k.','markersize',sizeofmark)
	plot([1 Tmax],[1 1],'y')
	xlabel('prediction time T')
	ylabel('NRMSE(T)')
	title([tittxt,' NRMSE(T) for prediction with AR(',int2str(m),...
            '), n=',int2str(n),' nlast=',int2str(nlast)])
end
