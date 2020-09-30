function xV = generateARMAts(phiV,thetaV,n,sdnoise)
% xV = generateARMAts(phiV,thetaV,n,sdnoise)
% Generate an ARMA(p,q) time series of length 'n' with Gaussian input noise.
% Note that phiV = [phi(0) phi(1) ... phi(p)]' and phi(0) is the constant
% term, and thetaV = [theta(1) ... theta(q)]'. 
% sdnoise is the SD of the input noise (if left out then sdnoise=1).
% The generating ARMA(p,q) process reads
% x(t) = phi(0) + phi(1)*x(t-1) + ... + phi(p)*x(t-p) + 
%        +z(t) - theta(1)*z(t-1) + ... - theta(q)*z(t-p), 
% z(t) ~ WN(0,sdnoise^2)

if nargin==2
    sdnoise=1;
end
if isempty(phiV)
    phiV = 0;
    p=0;
else
    p = length(phiV)-1;
end
if isempty(thetaV)
    q = 0;
else
    q = length(thetaV);
end
pq = max(p,q);
ntrans = 100+pq;

phiV = phiV(:);
thetaV = thetaV(:);
if p>0
    rootarV = roots([1;-phiV(2:end)]);
    if any(abs(rootarV)>=1)
        fprintf('The AR(%d) part of the process is not stationary.\n',p);
    end
end
if q>0
    rootmaV = roots([1;-thetaV(2:end)]);
    if any(abs(rootmaV)>=1)
        fprintf('The MA(%d) part of the process is not reversible.\n',q);
    end
end
x0V = sdnoise*randn(pq,1);
zV = randn(n+ntrans,1) * sdnoise;
xV = NaN*ones(n+ntrans,1);
xV(1:pq) = x0V;
if p==0
    for i=pq+1:n+ntrans
        xV(i)=phiV(1)+zV(i)-thetaV'*flipud(zV(i-q:i-1));
    end
elseif q==0    
    for i=pq+1:n+ntrans
        xV(i)=phiV(1)+phiV(2:p+1)'*flipud(xV(i-p:i-1))+zV(i);
    end
else
    for i=pq+1:n+ntrans
        xV(i)=phiV(1)+phiV(2:p+1)'*flipud(xV(i-p:i-1))+zV(i)-thetaV'*flipud(zV(i-q:i-1));
    end
end
xV = xV(ntrans+1:n+ntrans);
