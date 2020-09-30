function fnnM = falsenearest(xV,tau,mmax,escape,theiler,tittxt)
% fnnM = falsenearest(xV,tau,mmax,escape,theiler,tittxt)
% FNN computes the percentage of false nearest neighbors for a range of
% embedding dimensions starting from 1 to 'mmax' embedding dimensions. 
% For the search of nearest neighbors the Matlab functions for k-d-tree
% search are used.
% INPUT 
%  xV       : vector of the scalar time series
%  tau      : the delay time. If empty, tau=1
%  mmax     : the maximum embedding dimension.
%  escape   : A factor of escaping from the neighborhood. Default=10.
%  theiler  : the Theiler window to exclude time correlated points in the
%             search for neighboring points. Default=0.
%  tittxt   : string to be displayed in the title of the figure 
%             if not specified, no plot is made
% OUTPUT: 
%  fnnM     : A matrix of two columns and 'mmax' rows, where the embedding 
%             dimension is in the first column and the percentage of fnn in
%             the second column.
fthres = 0.1; % A factor of the data SD to be used as the maximum radius 
              % for searching for valid nearest neighbor.
propthres = 0.1; % Limit for the proportion of valid points, i.e. points 
                 % for which a nearest neighbor was found. If the proporion 
                 % of valid points is beyond this limit, do not compute
                 % FNN.              
thresh = 0.01;
n = length(xV);
if nargin==5
    tittxt = [];
elseif nargin==4
    tittxt = [];
    escape=10;
elseif nargin==3
    tittxt = [];
    escape=10;
    theiler=0;
end
if isempty(tau), tau=1; end
if isempty(escape), escape=10; end
if isempty(theiler) || theiler<0, theiler=0; end
theiler = round(theiler);
if n<=2*(theiler+1)
    error('Too large Theiler window=%d for the time series length n=%d',theiler,n); 
end
% Rescale to [0,1] and add infinitesimal noise to have distinct samples
xmin = min(xV);
xmax = max(xV);
xV = (xV - xmin) / (xmax-xmin);
xV = AddNoise(xV,10^(-10));
rthresmax = fthres*std(xV); % The maximum distance to look for nearest neighbor
fnncountV = NaN*ones(mmax,1);
m=1;
nextm = 1;
while m<=mmax && nextm
    nvec = n-m*tau; % to be able to add the component x(nvec+tau) for m+1 
    xM = NaN*ones(nvec,m);
    for i=1:m
        xM(:,m-i+1) = xV(1+(i-1)*tau:nvec+(i-1)*tau);
    end
    % k-d-tree data structure of the training set for the given m
    kdtreeS = KDTreeSearcher(xM);
    % For each target point, find the two nearest neighbors, the fitst one 
    % is itself since tha target point is included in the training set.
    if theiler == 0
        [neiindM,neidisM]=knnsearch(kdtreeS,xM,'K',2);
        idxV = neiindM(:,2); 
        distV = neidisM(:,2);
    else
        [neiindM,neidisM]=knnsearch(kdtreeS,xM,'K',2*theiler+2);
        idxV = NaN(nvec,1);
        distV = NaN(nvec,1);
        for ivec = 1:nvec
            i1 = max(ivec-theiler,1);
            i2 = min(ivec+theiler,nvec);
            ineiV = find(neiindM(ivec,:)<i1 | neiindM(ivec,:)>i2); 
            inei = ineiV(1);
            idxV(ivec) = neiindM(ivec,inei);
            distV(ivec) = neidisM(ivec,inei);
        end
    end
    iV = find(distV< rthresmax*sqrt(m));
    if isempty(iV)
        nextm = 0;
        nproper = 0;
    else
        nproper = length(iV);
    end
    % Compute fnn only if there is a sufficient number of target points 
    % having nearest neighbor (in R^m) withing the threshold distance
    if nproper>propthres*nvec
        nnfactorV = 1+(xV(iV+m*tau)-xV(idxV(iV)+m*tau)).^2./distV(iV).^2;
        fnncountV(m) = length(find(nnfactorV > escape^2))/nproper;
    end
    m = m+1;
end
fnnM = [[1:mmax];fnncountV']';
if ~isempty(tittxt)
	figno = gcf;
	figure(figno)
	clf
	plot([1:mmax]',fnncountV,'.-k')
    hold on
	plot([1 mmax],thresh*[1 1],'c--')
	xlabel('m')
	ylabel('FNN(m)')
	title([tittxt,' FNN (\tau=',int2str(tau),' w=',int2str(theiler),...
        ' f=',int2str(escape),'), n=',int2str(n)])
end
