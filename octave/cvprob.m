% convection-diffusion test problems
% from "On the convergence..." by Frommer and Szyld

% a coef of d^2_x, b for d^2_y, c for d_x, d for d_y
% or positive diff consts a,b
% velocity field c,d
function A = cdprob(a,b,v,mu,tau)



N = size(a,1);
h = 1/(N+1);
n = N^2;

if nargs < 5
  tau = 0.5*h*h;
end

% flatten fields
a = a(:);
b  = b(:);
v     = v(:);
mu    = mu(:);

% x-direction neighbors
ex = ones(N,1);
Tx = spdiags([-ex ex],[-1 1],N,N);
Ix = speye(N);

% y-direction neighbors
ey = ones(N,1);
Ty = spdiags([-ey ey],[-1 1],N,N);
Iy = speye(N);

% diffusion coefficients at half steps
a_p = zeros(n,1);
a_m = zeros(n,1);
b_p  = zeros(n,1);
b_m  = zeros(n,1);

% x-direction half points
a_p(1:end-1) = 0.5*(a(1:end-1) + a(2:end));
a_m(2:end)   = 0.5*(a(2:end)   + a(1:end-1));

% y-direction half points
b_p(1:end-N) = 0.5*(b(1:end-N) + b(1+N:end));
b_m(1+N:end) = 0.5*(b(1+N:end) + b(1:end-N));

% x-diffusion
Bx_diff = spdiags(a_p/h^2, 1, n, n) + ...
          spdiags(a_m/h^2,-1, n, n) - ...
          spdiags((a_p+a_m)/h^2, 0, n, n);

% y-diffusion
By_diff = spdiags(b_p/h^2,  N, n, n) + ...
          spdiags(b_m/h^2, -N, n, n) - ...
          spdiags((b_p+b_m)/h^2, 0, n, n);

% x-convection
Bx_conv = spdiags(v/(2*h), 1, n, n) - ...
          spdiags(v/(2*h),-1, n, n);

% y-convection
By_conv = spdiags(mu/(2*h),  N, n, n) - ...
          spdiags(mu/(2*h), -N, n, n);

% assemble full operator
B = Bx_diff + By_diff + Bx_conv + By_conv;
A = speye(n) + (tau/2)*B;

end


%!test
%! N = 100;
%! h = 1/(N+1);
%! tau = 0.5*h^2;
%! [x,y] = meshgrid((1:N)*h,(1:N)*h);
%! z = sin(pi*x).*sin(pi*y);
%! z = z(:);
%! alpha_fun = @(x,y) 1;
%! beta_fun  = @(x,y) 1;
%! v_fun     = @(x,y) 0;
%! mu_fun    = @(x,y) 0;
%! a= alpha_fun(x,y);
%! b= beta_fun(x,y);
%! v= v_fun(x,y);
%! mu= mu_fun(x,y);
%! A = cdprob(a,b,v,mu,tau);
%! AzNum = A*z
%! lambda_h = (8/h^2)*sin(pi*h/2)^2;
%! AzExact = (1 + tau/2*lambda_h)*z
%! rerr = norm(AzNum - AzExact) / norm(AzExact);
%! assert(rerr,e-13)
