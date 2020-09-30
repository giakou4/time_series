function flag = myversionOld
% flag = myversionOld
% Gives out a code number depending on the version of matlab running.
% OUTPUT
% flag = 0, if release is up to 2006b
% flag = 1, if release is from 2007a to 2010a
% flag = 2, if release is from 2010b onwards

% a -> 0
% b -> 0.5
rel1year = 2007.0;
rel2year = 2010.0;

relS = version('-release');
relyear = str2num(relS(1:4));
if relS(5)=='b'
    relyear = relyear+0.5;
end
if relyear < rel1year
    flag = 0;
elseif relyear < rel2year
    flag=1;
else
    flag=2;
end

