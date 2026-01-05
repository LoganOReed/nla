% gram schmidt
function [Q,R] = gs(A)
  [m,n] = size(A); % num cols
  Q = zeros(m,n);
  R = zeros(n,n);

  for j = 1:n
    qb = A(:,j);
    for i = 1:j-1
      R(i,j) = Q(:,i)' * A(:,j);
      qb = qb - R(i,j)*Q(:,i);
    end
    R(j,j) = norm(qb);
    Q(:,j) = qb / R(j,j);
  end
end

%!test
%! a = [0,-20,-14;3,27,-4;4,11,-2];
%! [q,r] = gs(a);
%! assert(size(q), size(a))
%! assert(size(r), [size(a,2), size(a,2)])

%!test
%! a = eye(10);
%! [q,r] = gs(a);
%! assert (eye(10), q*r);

%!test
%! a = [0,-20,-14;3,27,-4;4,11,-2];
%! [q,r] = gs(a);
%! assert(q*r,a,tol=1e-14)

