function plotrealpred(oriV,preM,legtxtM)
% plotrealpred(oriV,preM,legtxtM)
% PLOTREALPRED plots together predicted and actual sample points
% INPUTS:
%  oriV    : vector of actual scalar data 
%  preM    : matrix of predictions with different models (one model 
%           per column) 
%  legtxtM : a string matrix of the legends for each data set: 
%            first string for the original data set, second to last for
%            the data sets of predictions with different models at the 
%            clumn order in 'preM'.

sizeofmark = 6;
m = size(preM,2);
if nargin < 3
    legtxtM = 'real';
    for i=1:m
        legtxtM = str2mat(legtxtM,['predict-',int2str(i)]);
    end
end
symb1V = str2mat('''k''', '''r''', '''b''', '''g''', '''c''');
symb2V = str2mat('''k.''', '''r.''', '''b.''', '''g.''', '''c.''');
n = length(oriV);
iV = [1:n]';
figure(gcf)
clf
eval(['plot(iV,oriV,',symb1V(1,:),',''linewidth'',1)'])
if nargin==3
    legtxt = ['''real'','];
end
hold on
for i=1:m
    eval(['plot(iV,preM(:,i),',symb1V(i+1,:),',''linewidth'',1)'])
    legtxt = [legtxt,'''',legtxtM(i,:),''','];
end
eval(['plot(iV,oriV,',symb2V(1,:),')'])
for i=1:m
    eval(['plot(iV,preM(:,i),',symb2V(i+1,:),')'])
end
eval(['legend(',legtxt,'0)'])
xlabel('prediction time T')
ylabel('x(T)')
title('predictions of a time series')
