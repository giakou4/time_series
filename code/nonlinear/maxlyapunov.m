function [l1M,sdl1V] = maxlyapunov(xV,tau,mmin,mmax,tittxt,nitecal,itewidth)
% [l1V,sdl1V] = maxlyapunov(xV,tau,mmin,mmax,tittxt,nitecal,itewidth)
% This function computes the largest Lyapunov exponent for a given 
% time series 'xV', for a delay 'tau' and a range of embedding
% dimensions 'mmin,...,mmax'. The algorithm for the computation of the
% the largest Lyapunov exponent is 'lyap_k' in TISEAN and is presented
% in [H. Kantz, A robust method to estimate the maximal Lyapunov exponent 
% of a time series, Phys. Lett. A 185, 77 (1994)].
% INPUTS:
%  xV      : vector of a scalar time series
%  tau     : the delay time
%  mmin    : the minimum embedding dimension to compute the largest 
%            Lyapunov for. If mmin<2 then mmin=2.
%  mmax    : the maximum embedding dimension to compute the largest 
%            Lyapunov for, i.e. for embedding dimensions mmin,...,mmax
%  tittxt  : string to be displayed in the title
%  nitecal : number of first iterations to constrain the search for the  
%            characteristic slope that gives the maximum Lyapunov exponent.
%            Default is 20;
%  itewidth: sliding window of iterations to estimate local slopes and select
%            the one with the smallest fitting error. This is the maximum 
%            Lyapunov exponent. Default is 5;
% OUTPUT:
%  l1M     : matrix of size 'mmax-mmin+1'x2 of the largest Lyapunov exponents 
%            for the embedding dimensions m=mmin,...,mmax. These are computed
%            from the average of the estimated slopes for each epsilon (the
%            default number of 5 epsilon values in 'lyap_k' is adopted here).
%            The first column contains the average over epsilon values and the
%            second column the weighted average using the number of points found
%            to have enough neighbors.
%  sdl1V   : vector of length 'mmax-mmin+1' of the stndard deviation in the
%            estimation of the largest Lyapunov exponents for the embedding 
%            dimensions m=mmin,...,mmax. 
sizeofmark = 3;
if nargin==6
    itewidth = 5;
elseif nargin==5
    itewidth = 5;
    nitecal = 20;
elseif nargin==4
    itewidth = 5;
    nitecal = 20;
    tittxt = '';
end
if mmin<2
    mmin=2;
end
if mmax<2
    mmax=2;
end
neps = 5; % default number of epsilon values in lyap_k
nite = 50+1; % number of iterations in lyap_k
nminerr = round(nitecal/2);
mV = [mmin:mmax]';
nm = length(mV);
itecalV = [0:nitecal-1]';

symb1V = str2mat('''k-''', '''k--''', '''k-.''', '''k:''');
symb2V = str2mat('''b-''', '''b--''', '''b-.''', '''b:''');
symb3V = str2mat('''r-''', '''r--''', '''r-.''', '''r:''');
symbV = str2mat(symb1V, symb2V, symb3V);

save tmp.dat xV -ascii
eval(['!c:\tisean\lyap_k tmp.dat -d',int2str(tau),' -m',int2str(mmin),' -M',int2str(mmax),' -o tmp.lle'])

nneiM = NaN*ones(neps,nm);
jbestM = NaN*ones(neps,nm);
bbestM = NaN*ones(neps,nm);
nonemptyM = zeros(neps,nm);
for im=1:nm
    eval(['loglm',int2str(im),'M=NaN*ones(nitecal,neps);'])
end
epsV = NaN*ones(neps,1);
fid = fopen('tmp.lle','r');
for ieps=1:neps 
    % disp(['ieps =',int2str(ieps),'...'])
    for im=1:nm
        h1=fgets(fid);
        if im==1
            epsV(ieps) = sscanf(h1,'#epsilon= %f');
        end
        tmp = fgets(fid);
        if ~isempty(str2num(tmp(1)))
            tmp1 = sscanf(tmp,'%f %f %f');
            nonemptyM(ieps,im)=1;
            tmpM = fscanf(fid,'%f ',[3 nite])';
            tmpM = [tmp1'; tmpM];
            loglV = tmpM(1:nitecal,2);
            eval(['loglm',int2str(im),'M(:,ieps)=loglV;'])
            nneiM(ieps,im)=tmpM(1,3);
			sV = NaN*ones(nitecal-itewidth+1,1);
			for j=1:nitecal-itewidth+1
                [p,s]=polyfit(itecalV(j:j+itewidth-1),loglV(j:j+itewidth-1),1);
                bV(j)=p(1);
                sV(j)=s.normr;
			end
			[a,indV]=sort(sV); % The segments with smallest errorfit
			[bmax,imax] = max(bV(indV(1:nminerr)));
			[bmin,imin] = min(bV(indV(1:nminerr)));
			neibmaxV = []; % The indices of those with large slope
			neibminV = []; % The indices of those with small slope
			for k=1:nminerr
                if abs(bV(indV(k))-bmax)<abs(bV(indV(k))-bmin)
                    neibmaxV = [neibmaxV;indV(k)];
                else
                    neibminV = [neibminV;indV(k)];
                end
			end
			if length(neibmaxV)>1
                if bmax-bmin<2*std(bV(neibmaxV))
                    jbestM(ieps,im) = indV(1); % not really different slopes, choose the one with smallest errorfit
                else
                    jbestM(ieps,im) = neibmaxV(1); % choose the one with the smallest errorfit of the large slopes
                end
			elseif bmax-bmin<2*std(bV(neibminV))
                jbestM(ieps,im) = indV(1); % not really different slopes, choose the one with smallest errorfit
			else
                jbestM(ieps,im) = neibmaxV; % a single large slope 
			end
            bbestM(ieps,im)=bV(jbestM(ieps,im));
			% disp(['m=',int2str(mV(im)),'  jbest = ',int2str(jbestM(ieps,im)),...
            %        '  b(jbest) = ',num2str(bbestM(ieps,im))])
        end
	end
end
fclose(fid);
for im=1:nm
    figure(im)
    clf
    hold on
    legtxt = [];
    leastone = 'n';
    for ieps=1:neps
        if nonemptyM(ieps,im)==1
            eval(['plot(itecalV(jbestM(ieps,im):jbestM(ieps,im)+itewidth-1),loglm',...
                int2str(im),'M(jbestM(ieps,im):jbestM(ieps,im)+itewidth-1,ieps),',symbV(ieps,:),...
                ',''linewidth'',',int2str(sizeofmark),')'])
            legtxt = [legtxt,'''',int2str(ieps),' #x=',int2str(nneiM(ieps,im)),' \lambda_1=',num2str(bbestM(ieps,im),3),''','];
            leastone='y';
        end
    end    
    for ieps=1:neps
        eval(['plot(itecalV,loglm',int2str(im),'M(:,ieps),',symbV(ieps,:),')'])
    end    
    xlabel('iteration')
    ylabel('log(\lambda_1)')
    title([tittxt,' m=',int2str(mV(im))])
    if leastone=='y'
       % eval(['legend(',legtxt,'0)'])
    else
        ax = axis;
        text(ax(1)+0.3*(ax(2)-ax(1)),ax(3)+0.5*(ax(4)-ax(3)),['no calculations for m=',int2str(mV(im))])
    end
end
l1M = NaN*ones(nm,2);
l1M(:,1) = (nanmean(bbestM))';
sdl1V = (nanstd(bbestM))';
for im=1:nm
    if any(nonemptyM(:,im))
        l1M(im,2) = nansum(nneiM(:,im).*bbestM(:,im))/nansum(nneiM(:,im));
    end
end
figure(nm+1)
clf
errorbar(mV,l1M(:,1),sdl1V)
xlabel('m')
ylabel('\lambda_1(m)')
title([tittxt,' Lyapunov exponents'])
!del tmp.lle
