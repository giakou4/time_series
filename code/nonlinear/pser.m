function pser(xV, tittxt, nseg, type, tau_s, ignore)
% pser(xV, tittxt, nseg, type, tau_s, ignore)
% PSER gives the time history diagram of a time series in one
% figure, optionally containing subplots for successive segments.
% INPUTS:
%  xV      : vector of a scalar time series
%  tittxt  : string to be displayed in the title
%  nseg    : number of successive segments plotted on each other. 
%  type    : if 'd' (for discrete) then data points are displayed 
%            with dots, if 'c' lines are used and if 'b' lines and 
%            circles are used. 
%  tau_s   : sampling time in order to specify the real time on x-axis.
%            For dicrete-type data tau_s should be assigned to one. 
%  ignore  : if specified and it is 'y' (yes) then the last samples 
%            (equal to length(xV) mod nseg) are ignored. 

sizeofmark = 6;
if nargin < 6
    ignore = 'y';
end

n = length(xV);
n1 = floor(n / nseg);
if ignore ~= 'y' & n / nseg ~= n1  
    nsub = nseg+1;
else
    nsub = nseg;
end
xmin = min(xV);
xmax = max(xV);
dmax = xmax - xmin;
xmin = xmin - 0.05*dmax;
xmax = xmax + 0.05*dmax;
if type == 'd' & tau_s ~= 1
    error('The type is discrete and the sampling time different from 1.')
end
indV = [1:n]' * tau_s;
  
figure(gcf)
clf
for i=1:nseg
    if type == 'd'
        subplot(nsub,1,i), plot(indV((i-1)*n1+1:i*n1), xV((i-1)*n1+1:i*n1), 'k.')
    elseif type == 'c'
        subplot(nsub,1,i), plot(indV((i-1)*n1+1:i*n1), xV((i-1)*n1+1:i*n1), 'k')
    else
        subplot(nsub,1,i), plot(indV((i-1)*n1+1:i*n1), xV((i-1)*n1+1:i*n1), 'k')
        hold on
        subplot(nsub,1,i), plot(indV((i-1)*n1+1:i*n1), xV((i-1)*n1+1:i*n1), 'k.','markersize',sizeofmark)
    end 
    axis([indV((i-1)*n1+1) indV(i*n1) xmin xmax])
end
if nsub > nseg
    if type == 'd'
        subplot(nsub,1,nsub), plot(indV(nseg*n1+1:n), xV(nseg*n1+1:n), '.')
    elseif type == 'c'
        subplot(nsub,1,nsub), plot(indV(nseg*n1+1:n), xV(nseg*n1+1:n))
    else
        subplot(nsub,1,nsub), plot(indV(nseg*n1+1:n), xV(nseg*n1+1:n))
        subplot(nsub,1,nsub), plot(indV(nseg*n1+1:n), xV(nseg*n1+1:n), '.','markersize',sizeofmark)
    end  
    axis([indV(nseg*n1+1) (nseg+1)*n1*tau_s xmin xmax])
end
subplot(nsub,1,1), title([tittxt, ': time history, ', int2str(n), ' data'])
subplot(nsub,1,nsub), xlabel(['time index k (sampling time =', num2str(tau_s,2), ')'])
subplot(nsub,1,nsub), ylabel('x(k)') 