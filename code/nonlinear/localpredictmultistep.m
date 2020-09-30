function preV = localpredictmultistep(xV,n1,tau,m,Tmax,nnei,q,tittxt)
% preV = localpredictmultistep(xV,n1,tau,m,Tmax,nnei,q,tittxt);
% LOCALPREDICTMULTISTEP makes iterative multi-step ahead local predictions 
% using the zeroth order (average mapping or nearest neighbor mappings if 
% only one neighbor is chosen) or local linear approach. The k-d-tree data
% structure as implemented in Matlab is utilized to speed up computation 
% time in the search of neighboring points.
% INPUTS:
%  xV      : vector of the scalar time series
%  n1      : size of training set, number of samples of the segment of xV 
%            (x(1),x(2),...,x(n1)) from which the nearest neighbors are
%            formed to be used for local prediction. Note that for time n1 
%            the first target point is formed to make predictions for.
%  tau     : the delay time.
%  m       : the embedding dimension.
%  Tmax    : the prediction horizon, prediction are made for T=1...Tmax 
%            steps ahead.
%  nnei    : number of nearest neighbors to be used in the local model.
%  q       : the truncation parameter for a normalization of the local
%            linear model if specified (to project the parameter space of 
%            the model, using Principal Component Regression, PCR, locally).
%            if q>=m -> Ordinary Least Squares, OLS (standard local linear 
%            model, no projection)
%            if 0<q<m -> PCR(q)
%            if q=0 -> local average model (if in addition nnei=1 ->
%            then the zeroth order model is applied)
%  tittxt  : string to be displayed in the title of the figure (if omitted
%            or blank, no figure is generated).
% OUTPUT: 
%  preV    : vector of length Tmax of the predicted values 
%            x(n1+1),x(n1+2),...,x(n1+Tmax)
% The actual and predicted values are plotted (provided that the length of
% 'xV' is at least 'n1+Tmax').
sizeofmark = 6; 
if nargin==7
    tittxt = [];
elseif nargin==6
    tittxt = [];
    q=0;
elseif nargin==5
    tittxt = [];
    q=0;
    nnei=1;
elseif nargin==4
    tittxt = [];
    tau=1;    
    q=0;
    Tmax=1;
end
if isempty(tittxt), toplot=0; else toplot = 1; end
if isempty(tau), tau=1; end
if isempty(q), q=0; end
if isempty(nnei), nnei=1; end
if isempty(n1), n1=round(n/2); end
if isempty(Tmax), Tmax=1; end
if q>m, q=m; end
n = length(xV);
if n1>=n
    error('Size of the training set n1 must be smaller than the time series length.');
end
n1vec = n1-(m-1)*tau-1;
xM = NaN*ones(n1vec,m);
for i=1:m
    xM(:,m-i+1) = xV(1+(i-1)*tau:n1vec+(i-1)*tau);
end
kdtreeS = KDTreeSearcher(xM); % k-d-tree data structure of the training set
% For each target point, find neighbors, apply the linear models and keep
% track of the predicted values each prediction time.
preV = NaN*ones(Tmax,1);
winnowV = xV(n1-(m-1)*tau:n1)';
for T = 1:Tmax
    targV = winnowV(:,end:-tau:end-(m-1)*tau);
    neiindV=knnsearch(kdtreeS,targV,'K',nnei);
    neiM = xM(neiindV,:);
    yV = xV(neiindV+(m-1)*tau+1);
    if q==0 || nnei==1
        preV(T) = mean(yV);
    else
        mneiV = mean(neiM);
        my = mean(yV);
        zM = neiM - ones(nnei,1)*mneiV;
        [Ux, Sx, Vx] = svd(zM, 0);
        tmpM = Vx(:,1:q) * inv(Sx(1:q,1:q)) * Ux(:,1:q)';
        lsbV = tmpM * (yV - my);
        preV(T) = my + (targV-mneiV) * lsbV;
    end
    winnowV = [winnowV preV(T)];
end % for T
if toplot
    Tmax = min(length(preV),Tmax);
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
    if q == 0
        title([tittxt,' multi-step local average prediction, \tau=',int2str(tau),...
                ' m=',int2str(m),' n1=',int2str(n1),' K=',int2str(nnei)])
    else
        title([tittxt,' multi-step local linear prediction, \tau=',int2str(tau),...
                ' m=',int2str(m),' n1=',int2str(n1),' K=',int2str(nnei)])
    end
    legend('real','predicted',0)
end