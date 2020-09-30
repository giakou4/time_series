function [y] = AddNoise(x, lev);
% y = addnoise(x, lev);
% This function takes an input time series 'x' and adds Gaussian
% white noise to it with standard deviation 'lev' times the
% standard deviation of 'x'.

mn = mean(x);
sd = std(x);

sdno = sd * lev;

y = x + randn(length(x), 1) * sdno;
 
