function rhoV = ARcoef2autocorr(phiV,maxtau,display)
% rhoV = ARcoef2autocorr(phiV,maxtau,display)
% The function computes the autocorrelation function for lags up to a given
% lag 'maxtau' (or a default=5 if the second argument is not specified),
% for an autoregressive (AR) process determined in terms of the
% coefficients given in the array 'phiV' of length 'p'.
% X(t) = phiV(1) X(t-1) + ... +  phiV(p) X(t-p) + Z(t)
% The moments of the white noise innovation process Z(t) do not need to be 
% specified.
% The output 'rhoV' is an array of length 'maxtau' of the autocorrelation
% values for lags 1,2,...,maxtau.
% INPUTS 
% - phiV    : the p coefficients of the AR(p) process.
% - maxtau  : the maximum lag to compute autocorrelation for (if omitted,
%             then the default is maxtau=5.
% - display : if 1 then show a graph of the autocorrelation (if omitted no
%             figure will be generated).
% OUTPUTS
% - rhoV    : the array of size 'maxtau x 1' of the autocorrelation values.

if nargin==2
    display=0;
elseif nargin==1
    maxtau = 5;
    display = 0;
end
if isempty(maxtau)
    maxtau = 5;
end

p = length(phiV);
phiV = phiV(:);
aV = [1; -phiV];
gammaV = rlevinson(aV,1);
rhoV = gammaV(2:end)/gammaV(1);

for i=p+1:maxtau
    rhoV(i) = phiV' * rhoV(i-1:-1:i-p);
end
if display
    figure
    clf
    hold on
    for i=1:maxtau
        plot(i*[1 1],[0 rhoV(i)],'b','linewidth',2)
    end
    plot([0 maxtau+1],[0 0],'k')
    xlabel('\tau')
    ylabel('\rho(\tau)')
    artxt = ['X(t)='];
    for i=1:p
        if phiV(i) < 0
            artxt = [artxt,'-',num2str(abs(phiV(i)),2),'*X(t-',int2str(i),')'];
        elseif phiV(i) > 0
            artxt = [artxt,'+',num2str(phiV(i),2),'*X(t-',int2str(i),')'];
        end
    end
    artxt = sprintf('%s+Z(t)',artxt);
    title(sprintf('Autocorrelation of %s',artxt))
end



