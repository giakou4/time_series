function sV = movingaveragesmooth(xV,maorder) 
% sV = movingaveragesmooth(xV,maorder)
% MOVINGAVERAGESMOOTH makes a smoothing of the given time series 'xV' using
% a moving average filter of a given order 'maorder'. It slides the window
% of length 'maorder' over the whole time series, starting after the first 
% maorder/2 observations of the time series and ending before the last
% maorder/2 observations. So, the first and last maorder/2 components are 
% blank (NaN). If 'maorder' is an even number the two edges of the sliding
% window are weighted by 0.5.
% INPUTS 
% - xV      : vector of length 'n' of the time series
% - maorder : the maorder of the moving average filter
% OUTPUTS
% - sV      : vector of length 'n' of the smoothed time series

n = length(xV);
xV = xV(:);
sV = NaN*ones(n,1);
if maorder <=1
    return;
end
if mod(maorder,2)==0
    q = maorder/2;
    for i=q+1:n-q
        sV(i) = (0.5*xV(i-q)+sum(xV(i-q+1:i+q-1))+0.5*xV(i+q))/maorder;
    end
else
    q = (maorder-1)/2;
    for i=q+1:n-q
        sV(i) = sum(xV(i-q:i+q))/maorder;
    end
end    
