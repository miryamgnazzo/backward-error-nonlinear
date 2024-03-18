%
% Problem "beam" from NEP-PACK
%

n = 50000;

A0 = spdiags(ones(n, 1) * [1 -2 1], -1:1, n, n);
A0(end, end-1:end) = [ -n, n ];
A1 = sparse(n,n,1);
F = { speye(n), A0, A1 };
[V, L] = be_newton(F, @f, -linspace(0.1, 1, 4));
R = be_residual(F, @f, V, L);
norm(R, 'fro')

function [fv, fp] = f(lambda)
    fv = [ -lambda, 1.0,  exp(-lambda) ];
    fp = [ -1.0,    0.0, -exp(-lambda) ];
end