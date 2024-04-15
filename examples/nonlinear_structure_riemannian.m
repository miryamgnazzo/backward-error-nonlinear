% Quadratic eigenvalue problem, with possible nonlinear structures (low-rank) on the
% coefficients.
% It computes an approximation for the structured backward error, for 2
% approximate eigenpairs. It employs Riemannian optimization-based
% technique.
%
n = 100 %00;
k = 2;
l = 2;

 A0 = spdiags(ones(n, 1) * [1 -2 1], -1:1, n, n);
 U = randn(n, k);
 A1 = lowrank(U, -U);
 A2 = speye(n);
 
 F = { A0, A1, A2 };

% Approximation of 2 eigenpairs
[V, L] = be_newton(F, @f, [-1, -4]);

e = 1e-4;
dA0 = spdiags(randn(n, 3) * e, -1:1, n, n);
dA1 = struct; dA1.U = e * randn(n, k); dA1.V = -dA1.U;
dA2 = e * randn * speye(n);

A0t = A0 + dA0;
Ut = A1.U + dA1.U; A1t = lowrank(Ut, -Ut);
A2t = A2 + dA2;

Ft = { A0t, A1t, A2t };

% Norm of [A0 A1 A2] - [A0t A1t A2t]
Dtrue = cell(1, 3);
Dtrue{1} = A0 - A0t;
Dtrue{2} = lowrank([ A1.U, -A1t.U], [A1.V, A1t.V ]);
Dtrue{3} = A2 - A2t;

tic;
D = be_riemannian(Ft, @f, ...
    { 'sparse', 'low-rank', 'identity' }, V, L);
time = toc;

% Computed backward error
A0t2 = A0t + D{1}; 
A1t2 = lowrank([ A1t.U, D{2}.U ], [A1t.V, D{2}.V ]); A2t2 = A2t + D{3};

% Norm of the backward error
nrm = be_norm(D);
res = be_residual({ A0t2, A1t2, A2t2 }, @f, V, L);

fprintf('Backward error: %e, Residual: %e, Perturbation norm: %e\n', ...
    nrm, norm(res, 'fro'), be_norm(Dtrue));

function [fv, fp] = f(l)
    fv = [1.0, l, l.^2];
    fp = [0.0, 1.0, 2.0*l];
end