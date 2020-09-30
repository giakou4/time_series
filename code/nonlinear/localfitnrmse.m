function [nrmseV,preM] = localfitnrmse(xV,tau,m,Tmax,nnei,q,tittxt)
% [nrmseV,preM] = localfitnrmse(xV,tau,m,Tmax,nnei,q,tittxt)
% LOCALFITNRMSE makes fitting using a local model of zeroth order (average 
% mapping or nearest neighbor mappings if only one neighbor is chosen) or a 
% local linear model and computes the fitting error for T-step ahead. For 
% the search for neighboring points it uses the Matlab k-d-tree search.
% The fitting here means that predictions are made for all the points in
% the data set (in-sample prediction). The prediction error statistic 
% (NRMSE measure) for the T-step ahead predictions is the goodness-of-fit 
% statistic. 
% The state space reconstruction is done with the method of delays having 
% as parameters the embedding dimension 'm' and the delay time 'tau'. 
% The local prediction model is one of the following:
% Ordinary Least Squares, OLS (standard local linear model): if the 
% truncation parameter q >= m
% Principal Component Regression, PCR, project the parameter space of the 
% model to only q of the m principal axes: if 0<q<m
% Local Average Mapping, LAM: if q=0.
% The local region is determined by the number of neighbours 'nnei'. 
% The k-d-tree data structure is utilized to speed up computation time in 
% the search of neighboring points and the implementation of Matlab is 
% used. 
% INPUTS:
%  xV      : vector of the scalar time series
%  tau     : the delay time (usually set to 1).
%  m       : the embedding dimension.
%  Tmax    : the prediction horizon, the fit is made for T=1...Tmax steps
%            ahead.
%  nnei    : number of nearest neighbors to be used in the local model. 
%            If k=1,the nearest neighbor mapping is the fitted value. 
%            If k>1, the model as defined by the input patameter 'q' is
%            used. 
%  q       : the truncation parameter for a normalization of the local
%            linear model if specified (to project the parameter space of
%            the model, using Principal Component Regression, PCR, locally).
%            if q>=m -> Ordinary Least Squares, OLS (standard local linear
%                       model, no projection)
%            if 0<q<m -> PCR(q)
%            if q=0 -> local average model (if in addition nnei=1 ->
%            then the zeroth order model is applied)
%  tittxt  : string to be displayed in the title of the figure 
%            if not specified, no plot is made
% OUTPUT: 
%  nrmseV  : vector of length Tmax, the nrmse of the fit for T-mappings, 
%            T=1...Tmax.
%  preM    : the matrix of size nvec x (1+Tmax) having the fit (in-sample
%            predictions) for T=1,...,Tmax, for each of the nvec 
%            reconstructed points from the whole time series. The first
%            column has the time of the target point and the rest Tmax
%            columns the fits for T=1,...,Tmax time steps ahead.
sizeofmark = 10; 
n = length(xV);
if nargin==6
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
if isempty(tau), tau=1; end
if isempty(q), q=0; end
if isempty(nnei), nnei=1; end
if isempty(Tmax), Tmax=1; end
if q>m, q=m; end
if n<2*((m-1)*tau-Tmax)
    error('the length of the time series is too small for this data size');
end
nvec = n-(m-1)*tau-Tmax;
xM = NaN*ones(nvec,m);
for i=1:m
    xM(:,m-i+1) = xV(1+(i-1)*tau:nvec+(i-1)*tau);
end
kdtreeS = KDTreeSearcher(xM); % k-d-tree data structure of the training set
% For each target point, find neighbors, apply the local state space model
% and store the predicted values for each prediciton time T.
% For T=1, in-sample prediction, the target point is in the training set.
preM = NaN*ones(nvec,Tmax);
neiindM=knnsearch(kdtreeS,xM,'K',nnei+1);
neiindM(:,1) = [];
for i=1:nvec
    neiM = xM(neiindM(i,:),:);
    yV = xV(neiindM(i,:)+(m-1)*tau+1);
    if q==0 || nnei==1
        preM(i,1) = mean(yV);
    else
        mneiV = mean(neiM);
        my = mean(yV);
        zM = neiM - ones(nnei,1)*mneiV;
        [Ux, Sx, Vx] = svd(zM, 0);
        tmpM = Vx(:,1:q) * inv(Sx(1:q,1:q)) * Ux(:,1:q)';
        lsbV = tmpM * (yV - my);
        preM(i,1) = my + (xM(i,:)-mneiV) * lsbV;
    end
end    
% For T>1, the target point is not in the training set.
% For each target point, find neighbors, apply the linear models and keep
% track of the predicted values for each model and each prediction time.
if Tmax>1
    winnowM = NaN*ones(nvec,(m-1)*tau+1);
    for i=1:(m-1)*tau+1
        winnowM(:,i) = xV(1+(i-1):nvec+(i-1));
    end
    for T = 2:Tmax
        winnowM = [winnowM preM(:,T-1)];
        targM = winnowM(:,end:-tau:end-(m-1)*tau);
        neiindM=knnsearch(kdtreeS,targM,'K',nnei);
        for i=1:nvec
            neiM = xM(neiindM(i,:),:);
            yV = xV(neiindM(i,:)+(m-1)*tau+1);
            if q==0 || nnei==1
                preM(i,T) = mean(yV);
            else
                mneiV = mean(neiM);
                my = mean(yV);
                zM = neiM - ones(nnei,1)*mneiV;
                [Ux, Sx, Vx] = svd(zM, 0);
                tmpM = Vx(:,1:q) * inv(Sx(1:q,1:q)) * Ux(:,1:q)';
                lsbV = tmpM * (yV - my);
                preM(i,T) = my + (targM(i,:)-mneiV) * lsbV;
            end
        end  
    end % for T
end % if Tmax>1        
preM = [[1:nvec]+(m-1)*tau; preM']'; %Add the target point index before the 
                         % iterative predictions nrmseV = NaN*ones(Tmax,1);
nrmseV = NaN*ones(Tmax,1);
for T=1:Tmax
    nrmseV(T) = nrmse(xV(preM(:,1)+T),preM(:,T+1));
end
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
	title([tittxt,' NRMSE(T), fit LP(m=',int2str(m),...
            ' K=',int2str(nnei),' q=',int2str(q),'), n=',int2str(n)])
end
