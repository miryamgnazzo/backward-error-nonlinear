% Problem "beam" from NEP-PACK
% It employs a Riemannian optimization based technique,
% in order to provide an approximation of the structured backward error.
%
n = 100000;

A0 = speye(n);
A1 = spdiags(ones(n, 1) * [1 -2 1], -1:1, n, n);
A1(end, end-1:end) = [ -n, n ];
A2 = sparse(n,n,1);
F = { A0, A1, A2 };
[V, L] = be_newton(F, @f, -linspace(0.1, 1, 3));

e = 1e-4;
dA0 = e * randn * speye(n);
dA1 = spdiags(randn(n, 3) * e, -1:1, n, n);
dA2 = sparse(n,n, e * randn);

A0t = A0 + dA0;
A1t = A1 + dA1;
A2t = A2 + dA2;

Ft = { A0t, A1t, A2t };

% Norm of [A0 A1 A2] - [A0t A1t A2t]
Dtrue = cell(1, 3);
Dtrue{1} = A0 - A0t;
Dtrue{2} = A1 - A1t;
Dtrue{3} = A2 - A2t;

tic;
D = be_riemannian(Ft, @f, ...
    { 'identity', 'sparse', 'sparse' }, V, L);
time  = toc;

% Computed backward error
A0t2 = A0t + D{1}; 
A1t2 = A1t + D{2};
A2t2 = A2t + D{3};

% Norm of the backward error
nrm = be_norm(D);
res = be_residual({ A0t2, A1t2, A2t2 }, @f, V, L);

fprintf('Backward error: %e, Residual: %e, Perturbation norm: %e\n', ...
    nrm, norm(res, 'fro'), be_norm(Dtrue));

function [fv, fp] = f(lambda)
    fv = [ -lambda, 1.0,  exp(-lambda) ];
    fp = [ -1.0,    0.0, -exp(-lambda) ];
end