% CGW algo2 from Andrew,Szyld paper
function [x,R,iter,t] = algo3(A,b,tol,maxit)
%solves (M-N)x=b according to Widlund, A LANCZOS METHOD FOR A CLASS OF NONSYMMETRIC
%SYSTEMS OF LINEAR EQUATIONS

%x2 contains approx. solution at the last step
%R contains the norm of errors in each step
%t is the running time of recurrence
%iter : number of iterations

%Stops when tolerance or maxit is reached

%default tolerance
if nargin<3
    tol=1e-8;
end
if nargin<4
    maxit=100;
end

H = 0.5 * (A + A'); % Hermitian
S = 0.5 * (A - A'); % Skew-Hermitian
K = H\S;
[n,~]=size(A);


% NOTE: only store 3 for 3-term recurrence
vprev=zeros(n,1);
v = vprev;

xprev=zeros(n,1);
x=xprev;


aprev = 0;
a = 1;
pprev = 1;
p = 1;

rprev = b;
r = rprev;
R(1) = norm(r);
tic;
k = 1;
while (R(k)>tol) && (k<maxit)
  v = H \ r;
  p = r'*v;
  a = 1 + (p / pprev) * aprev;
  w = 1 / a;
  xnext = (1 - w)*xprev + w*x + w*v;
  rnext = (1 - w)*rprev - w*S*v;
  R(k+1) = norm(rnext);

  % Update the xs and vs so they line up with k
  aprev = a;
  pprev = p;

  rprev = r;
  r = rnext;

  xprev = x;
  x = xnext;

  vprev = v;

  k = k+1;
end
iter = k
t = toc;

end

%
% %! alpha_fun = @(x,y) 1;
% %! beta_fun  = @(x,y) 1;
% %! v_fun     = @(x,y) 0;
% %! mu_fun    = @(x,y) 0;
% %! a= alpha_fun(x,y);
% %! b= beta_fun(x,y);
% %! v= v_fun(x,y);
% %! mu= mu_fun(x,y);

%!test
%! N = 100;
%! h = 1/(N+1);
%! tau = 0.5*h^2;
%! [x,y] = meshgrid((1:N)*h,(1:N)*h);
%! z = x.*y.*(1-x).*(1-y);
%! z = z(:);
%! a = ones(N);
%! b = ones(N);
%! v = zeros(N);
%! mu = zeros(N);
%! A = cdprob(a,b,v,mu,tau);
%! [x,~,~,~] = algo1(A,A*z);
%! x=x
%! z=z
%! rerr = norm(z - x) / norm(x);
%! assert(rerr,0,tol=3e-3);
