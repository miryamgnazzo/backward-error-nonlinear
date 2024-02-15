
n = 32;
rng(123)

M = spdiags(ones(n, 1) * [1 2 1], -1:1, n, n);
K = spdiags(ones(n, 1) * [1 -2 1], -1:1, n, n);
D = spdiags(ones(n, 1) * [1 -1], 0:1, n, n);

epsilon = 1e-4;
dM = epsilon * spdiags(randn(n, 3), -1:1, n, n);
dK = epsilon * spdiags(randn(n, 3), -1:1, n, n);
dD = epsilon * spdiags(randn(n, 2), 0:1, n, n);

[V, L] = polyeig(K + dK, D + dD, M + dM);

s = 5;
V = V(:, 1:s);
L = L(1:s);

[D0, D1, D2] = sparse_QEP(K, D, M, V, diag(L));
% [D0, D1, D2] = struct_err_sparse(K, D, M, V, diag(L));

resnrm = norm(D2 * V * diag(L)^2 + D1 * V * diag(L) + D0 * V, 'fro');
fprintf('Residual norm for the backward error: %e\n', resnrm)

be_nrm = norm([ D0 D1 D2 ] - [ K D M ], 'fro');
fprintf('Structured backward error: %e\n', be_nrm);

fprintf('|| [ dM dD dK ] ||_F = %e\n', norm([ dM dD dK ], 'fro'));



