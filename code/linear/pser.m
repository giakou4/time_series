function pser(xV, nseg, tittxt, type, tau_s, ignore)
% pser(xV, tittxt, nseg, type, tau_s, ignore)
% PSER gives the time history diagram of a time series in one
% figure, optionally containing subplots for successive segments.
% INPUTS:
%  xV      : vector of a scalar time series
%  nseg    : number of successive segments plotted on each other. 
%  tittxt  : string to be displayed in the title
%  type    : if 'd' (for discrete) then data points are displayed 
%            with dots, if 'c' lines are used, and otherwise lines and 
%            dots are used. 
%  tau_s   : sampling time in order to specify the real time on x-axis.
%            For dicrete-type data tau_s should be assigned to one. 
%  ignore  : if specified and it is 1 (meaning yes) then the last samples 
%            (equal to length(xV) mod nseg) are ignored. 

if nargin == 5
    ignore = 1;
elseif nargin == 4
    ignore = 1;
    tau_s = 1;
elseif nargin == 3
    ignore = 1;
    tau_s = 1;
    type = 'b';
elseif nargin == 2
    ignore = 1;
    tau_s = 1;
    type = 'b';
    tittxt = [];
elseif nargin == 1
    ignore = 1;
    tau_s = 1;
    type = 'b';
    tittxt = [];
    nseg = 1;
end
if isempty(ignore)
    ignore = 1;
end
if isempty(tau_s)
    tau_s = 1;
end
if isempty(type)
    type = 'b';
end
if isempty(nseg)
    nseg = 1;
end

n = length(xV);
n1 = floor(n / nseg);
if ~ignore && n / nseg ~= n1  
    nsub = nseg+1;
else
    nsub = nseg;
end
xmin = min(xV);
xmax = max(xV);
dmax = xmax - xmin;
xmin = xmin - 0.05*dmax;
xmax = xmax + 0.05*dmax;
indV = [1:n]' * tau_s;
  
figure(gcf)
clf
for i=1:nseg
    if type == 'd'
        subplot(nsub,1,i), plot(indV((i-1)*n1+1:i*n1), xV((i-1)*n1+1:i*n1), 'k.')
    elseif type == 'c'
        subplot(nsub,1,i), plot(indV((i-1)*n1+1:i*n1), xV((i-1)*n1+1:i*n1), 'k')
    else
        subplot(nsub,1,i), plot(indV((i-1)*n1+1:i*n1), xV((i-1)*n1+1:i*n1), '.-k')
    end 
    axis([indV((i-1)*n1+1) indV(i*n1) xmin xmax])
end
if nsub > nseg
    if type == 'd'
        subplot(nsub,1,nsub), plot(indV(nseg*n1+1:n), xV(nseg*n1+1:n), 'k.')
    elseif type == 'c'
        subplot(nsub,1,nsub), plot(indV(nseg*n1+1:n), xV(nseg*n1+1:n), 'k')
    else
        subplot(nsub,1,nsub), plot(indV(nseg*n1+1:n), xV(nseg*n1+1:n),'.-k')
    end  
    axis([indV(nseg*n1+1) (nseg+1)*n1*tau_s xmin xmax])
end
subplot(nsub,1,1), title([tittxt, ': time history, ', int2str(n), ' data'])
subplot(nsub,1,nsub), xlabel(['t (sampling time =', num2str(tau_s,2), ')'])
subplot(nsub,1,nsub), ylabel('x(t)') 