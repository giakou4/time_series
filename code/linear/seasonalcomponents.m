function sV = seasonalcomponents(xV,per) 
% sV = seasonalcomponents(xV,per)
% SEASONALCOMPONENTS computes the periodic time series comprised of repetetive
% patterns of seasonal components given a time series and the season
% (period).
% INPUTS 
% - xV      : vector of length 'n' of the time series
% - per     : the season (period)
% OUTPUTS
% - sV      : vector of length 'n' of the time series of seasonal components.

n = length(xV);
xV = xV(:);
monV = NaN*ones(per,1);
sV = NaN*ones(n,1);
for j=1:per
    monV(j) = mean(xV(j:per:n));
end
monV = monV - mean(monV);
for j=1:per
    sV(j:per:n) = monV(j)*ones(length([j:per:n]),1);
end
