%
% Hadeler problem from the NLEVP collection
%

[F, f,~,~,SOL] = nlevp('fiber');

% The first two eigenvalues are good initial estimate for n = 8
% [V, L] = be_newton(F, f, SOL.eval);
[V, L] = be_newton(F, f, 7e-7);
%keyboard
%V=SOL.evec;
%L=SOL.eval;

% Perturb the problem
Ft = be_perturb(F, 1e-3);

% Compute the backward error, and compare with the estimates
D = be_unstructured(Ft, f, V, L);

fprintf('Norm of the Backward error: %e\n', be_norm(D));

for j = 1 : 3
    bnd = be_unstructured_bound(j, Ft, f, V, L);
    fprintf('BNDTYPE = %d, %e <= %e\n', j, be_norm(D), bnd);
end
