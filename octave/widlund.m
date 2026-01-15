function [x2,R,iter,t] = widlund(M,N,b,tol,maxit,L)
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

A=M-N;
[n,~]=size(A);


x1=zeros(n,1);%taking the initial vectors as zero
x2=zeros(n,1);
w=1;
k=2;
r=b;
nb=norm(b);

R(1)=1;% needs to be modified if the initial vectors are not zero
R(2)=1;

a(2)=1;

tic;
while (R(k)>tol) && (k<maxit)
    
    
    V_k=L'\(L\r);
    p(k)=V_k'*(M*V_k);
    V_k=V_k/sqrt(p(k));
    
    if k>2
        a(k)=(1+(p(k)/p(k-1))*a(k-1));
        w=1/a(k);
    end
    
    
    temp=x1+w*(sqrt(p(k))*V_k+x2-x1);
    
    r=b-A*temp;
    
    R(k+1)=norm(r)/nb;
    
    k=k+1;
    x1=x2;
    x2=temp;
end
iter=k-1;
t=toc;
R=R(2:end);
end

