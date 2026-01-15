% testing stuff

N = 100;
h = 1/(N+1);
tau = 0.5*h^2;
[x,y] = meshgrid((1:N)*h,(1:N)*h);
z = x.*y.*(1-x).*(1-y);
z = z(:);
a = ones(N) ./ 2;
b = ones(N) ./2;
v = zeros(N);
mu = zeros(N);
A = cdprob(a,b,v,mu,tau);

[x,R,it,t] = algo4(A,A*z,1e-8,200);
iterations = it
time = t
rerr = norm(z - x) / norm(z)
lastResidual = R(end)

% Create a figure
figure;
%
% Plot residuals vs iterations
semilogy(1:it, R, '-o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Iteration');
ylabel('Residual');
title('Algo4: Residuals over Iterations');
grid on;
%
% Save the figure as a PNG
% print('algo1_simple.png', '-dpng', '-r300'); % -r300 sets resolution to 300 dpi

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
