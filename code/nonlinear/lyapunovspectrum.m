function [lM,KYV] = lyapunovspectrum(xV,tau,mmax,k,tittxt)
% [lM,KYV] = lyapunovspectrum(xV,tau,mmax,k,tittxt)
% This function computes the Lyapunov spectrum for a given 
% time series 'xV', for a delay 'tau' and a range of embedding
% dimensions '2,...,mmax'. 
% INPUTS:
%  xV      : vector of a scalar time series
%  tau     : the delay time
%  mmax    : the maximum embedding dimension to compute Lyapunov 
%            exponents for, i.e. for embedding dimensions 2,...,mmax
%  k       : the number of neighbors to use for each local map
%  tittxt  : string to be displayed in the title
% OUTPUT:
%  lM      : matrix of size (mmax-1)*(mmax-1), where each row i has the Lyapunov
%            exponents for the embedding dimension m=i+1 (because it starts at m=2).
%  KYdim   : the vector of length (mmax-1) of the estimated Kaplan-Yorke dimensions.
%            The Kaplan-Yorke dimension is useful (and then it is a good estimate 
%            of the fractal dimension) only when the embedding dimension is 
%            proper (i.e. the integer larger than the fractal dimension).
sizeofmark = 6;
if nargin<4
    tittxt = '';
end

n = length(xV);
if n > 500
    k = round(n/100);
else
    k=5;
end
figure(gcf)
clf
for im = 1:mmax-1
    m = im+1;
    disp(['Running for m=',int2str(m),' ...'])
    [xM] = embeddelays(xV, m, tau);
    save tmp.dat xM -ascii
    eval(['!c:\tisean\lyap_spec tmp.dat -I -m',int2str(m),',1 -k',int2str(k),' -V0 -o tmp.out'])
    fid = fopen('tmp.out', 'r');
    tmp1S = fgets(fid);
    while tmp1S(1)~='#'
        tmpS = tmp1S;
        tmp1S = fgets(fid);
    end
    tmpV = sscanf(tmpS,'%f ', m+1);
    lM(im,1:m) = tmpV(2:m+1)';
    fgets(fid);
    fgets(fid);
    fgets(fid);
    KYV(im) = fscanf(fid,'#estimated KY-Dimension= %f');
    fclose(fid);
    plot([1:m],lM(im,1:m))
    hold on
    plot([1:m],lM(im,1:m),'.','markersize',sizeofmark)
    !del tmp.out
end
!del tmp.dat
