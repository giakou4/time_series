function sV = movingaverageseasonal(xV,maorder) 
% sV = movingaverageseasonal(xV,maorder)
% MOVINGAVERAGESEASONAL computes the periodic time series comprised of
% repetetive patterns of seasonal components given a time series and the
% season (period). First a moving average smoothing is applied, and then
% this is subtracted from the time series. The moving average smoothing is
% computing by the average of the values in a time window of length
% 'maorder' sliding across the whole time series, starting after the first 
% maorder/2 observations of the time series and ending before the last
% maorder/2 observations. So, the first and last maorder/2 components are 
% blank (NaN). If 'maorder' is an even number the two edges of the sliding
% window are weighted by 0.5.
% INPUTS 
% - xV      : vector of length 'n' of the time series
% - maorder : the season or period, which is also the order of the moving
%             average filter. 
% OUTPUTS
% - sV      : vector of length 'n' of the time series of seasonal
%             components. If maorder <=1, then sV has n zeros.

n = length(xV);
xV = xV(:);
wV = NaN*ones(n,1);
if maorder <=1
    sV = zeros(n,1);
    return;
end
if mod(maorder,2)==0
    q = maorder/2;
    for i=q+1:n-q
        wV(i) = (0.5*xV(i-q)+sum(xV(i-q+1:i+q-1))+0.5*xV(i+q))/maorder;
    end
else
    q = (maorder-1)/2;
    for i=q+1:n-q
        wV(i) = sum(xV(i-q:i+q))/maorder;
    end
end
seasonV = NaN*ones(maorder,1);
for i=1:maorder
    seasonV(i) = mean(xV(q+i:maorder:n-q) - wV(q+i:maorder:n-q));
end
seasonV = seasonV - mean(seasonV);
sV = NaN*ones(n,1);
for i=1:maorder
    sV(q+i:maorder:n-q) = seasonV(i)*ones(length([q+i:maorder:n-q]),1);
end

