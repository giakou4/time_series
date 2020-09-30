function [preV] = linearpredictmultistep(xV,n1,m,Tmax,tittxt);
% [preV] = linearpredictmultistep(xV,n1,m,Tmax,tittxt);
% LINEARPREDICTMULTISTEP makes multi-step ahead predictions 
% using an AR model.
% INPUTS:
%  xV      : vector of the scalar time series
%  n1      : size of training set, number of samples of the segment of xV 
%            (x(1),x(2),...,x(n1)) to which the AR model is fitted.
%  m       : the order of the AR model
%  Tmax    : the prediction horizon, prediction are made for T=1...Tmax steps
%            ahead.
%  tittxt  : string to be displayed in the title of the figure (if plotted)
% OUTPUT: 
%  preV    : vector of length Tmax of the predicted values 
%            x(n1+1),x(n1+2),...,x(n1+Tmax)
% The actual and predicted values are plotted (provided that the length of
% 'xV' is at least 'n1+Tmax').
sizeofmark = 6; 

if nargin==4
    tittxt = [];
end
mx = mean(xV(1:n1));
yV = xV(1:n1)-mx;
th = ar(yV,m);
aV = get(th,'a');
aV = -aV(2:m+1);
preV = NaN*ones(m+Tmax,1);
preV(1:m)=xV(n1-m+1:n1)-mx;
for T=m+1:m+Tmax
    preV(T)=aV*preV(T-1:-1:T-m);
end
preV = preV(m+1:m+Tmax);
preV = preV+mx;
iV = [n1+1:n1+Tmax]';
if length(xV)<n1+Tmax
    i2V = [n1+1:length(xV)]';
    oriV = xV(i2V);
elseif length(xV)==n1
    i2V = [];
    oriV = [];
else
    i2V = iV;
    oriV = xV(i2V);
end
if ~isempty(tittxt)
    figno = gcf;
    figure(figno)
    clf
    plot(i2V,oriV,'k')
    hold on
    plot(iV,preV,'r')
    plot(i2V,oriV,'k.','markersize',sizeofmark)
    plot(iV,preV,'r.','markersize',sizeofmark)
    xlabel('T')
    ylabel('x(t+T)')
    title([tittxt,' multi-step linear prediction, m=',int2str(m),' n1=',int2str(n1)])
    legend('real','predicted',0)
end