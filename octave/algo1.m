% CGW algo1 from Andrew,Szyld paper
function [x2,R,iter,t] = algo1(M,N,b,tol,maxit,L)
%solves (M-N)x=b according to Widlund, A LANCZOS METHOD FOR A CLASS OF NONSYMMETRIC
%SYSTEMS OF LINEAR EQUATIONS

%x2 contains approx. solution at the last step
%R contains the norm of errors in each step
%t is the running time of recurrence
%iter : number of iterations

%Stops when tolerance or maxit is reached

%default tolerance
if nargin<4
    tol=1e-8;
end
if nargin<5
    maxit=100;
end

H = 0.5 * (A + A'); % Hermitian
S = 0.5 * (A - A'); % Skew-Hermitian
K = H\S;
[n,~]=size(A);


% NOTE: only store 3 for 3-term recurrence
v1 = H \ b;
v2=zeros(n,1);
v3=zeros(n,1);

x1=zeros(n,1);%taking the initial vectors as zero
x2=v1;
x3=zeros(n,1);

% NOTE: Stored as list for visualization
r(1) = b;


% a(1) = 1;
b(1) = 0;


tic;
% NOTE: for v/x, 1 is prev, 2 is current, 3 is next
% so 1 = k-1, 2 = k, 3=k+1
k = 1;
while (r(k)>tol) && (k<maxit)
  p(k) = v2'*H*v2;
  if k == 1
    a(k) = 1;
  else
    a(k) = 1 + (p(k) / p(k-1)) * a(k-1)
  end
  w = 1 / a(k);
  if k == 1
    x3 = v1; 
  else
    x3 = (1 - w)*x1 + w*x2 + w*v2;
  end
  r(k+1) = b - A*x3;
  if k == 1
    v2 = v1;
    v3 = w*K*v2;
  else
    v3 = w*K*v2 + (1 - w)*v1
  end

  % Update the xs and vs so they line up with k
  x1 = x2;
  x2 = x3;

  v1 = v2;
  v2 = v3;

  k = k+1
end
t = toc;

end

