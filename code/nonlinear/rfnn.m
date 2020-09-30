n = 1000;
ts = 0.01;
tau = 5;
mmax = 10;
escape = 10;
theiler = 0;

xM = lorenzxyz(n,ts);
xV = xM(:,1);
xV = AddNoise(xV,0.1);
tic;
[fnn1M,mdist1V,sddist1V] = falsenearestM(xV,tau,mmax,escape,theiler,'Matlab');
t1 = toc;
tic;
[fnn2M,mdist2V,sddist2V] = falsenearest(xV,tau,mmax,escape,theiler,'kdtree');
t2 = toc;
fprintf('Time for Matlab=%f, for kdree=%f \n',t1,t2);
[fnn1M fnn2M]
