function rhoV = MAcoef2autocorr(thetaV,display)
% rhoV = MAcoef2autocorr(thetaV,maxtau,display)
% The function computes the autocorrelation function for a moving average
% (MA) process determined in terms of the coefficients given in the array 
% 'thetaV' of length 'q'. 
% X(t) = Z(t) - thetaV(1) Z(t-1) - ... -  thetaV(q) X(t-q) 
% The moments of the white noise innovation process Z(t) do not need to be 
% specified. 
% The output 'rhoV' is an array of length 'q' of the autocorrelation
% values for lags 1,2,...,q (for larger lags the autocorrelation is 0).
% INPUTS 
% - thetaV  : the q coefficients of the MA(q) process.
% - display : if 1 then show a graph of the autocorrelation (if omitted no
%             figure will be generated).
% OUTPUTS
% - rhoV    : the array of size 'q x 1' of the autocorrelation values.

if nargin==1
    display = 0;
end

q = length(thetaV);
thetaV = thetaV(:);
rhoV = NaN*ones(q,1);
denom = 1+sum(thetaV.^2);
for i=1:q
    numer = -thetaV(i);
    if q-i>0
        numer = numer+sum(thetaV(1:q-i).*thetaV(i+1:q));
    end
    rhoV(i)=numer/denom; 
end
if display
    figure
    clf
    hold on
    for i=1:q
        plot(i*[1 1],[0 rhoV(i)],'b','linewidth',2)
    end
    plot([0 q+1],[0 0],'k')
    xlabel('\tau')
    ylabel('\rho(\tau)')
    matxt = ['X(t)=Z(t)'];
    for i=1:q
        if thetaV(i) < 0
            matxt = [matxt,'+',num2str(abs(thetaV(i)),2),'*Z(t-',int2str(i),')'];
        elseif thetaV(i) > 0
            matxt = [matxt,'-',num2str(thetaV(i),2),'*Z(t-',int2str(i),')'];
        end
    end
    title(sprintf('Autocorrelation of %s',matxt))
end



