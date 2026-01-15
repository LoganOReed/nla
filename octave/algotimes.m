% testing stuff

N = 100;
h = 1/(N+1);
tau = 0.5*h^2;
[x,y] = meshgrid((1:N)*h,(1:N)*h);
z = x.*y.*(1-x).*(1-y);
z = z(:);
a = ones(N) ./ 2;
b = ones(N) ./ 2;
v = zeros(N);
mu = zeros(N);
A = cdprob(a,b,v,mu,tau);
[~,~,~,t1] = algo1(A,A*z,1e-8,400);
[~,~,~,t2] = algo2(A,A*z,1e-8,400);
[~,~,~,t3] = algo3(A,A*z,1e-8,400);
[~,~,~,t4] = algo4(A,A*z,1e-8,400);
[~,~,~,t5] = algo5(A,A*z,1e-8,400);
[~,~,~,t6] = algo6(A,A*z,1e-8,400);

tFourOneDiff = t1-t4
tFiveTwoDiff = t2-t5
tSixThreeDiff = t3-t6
