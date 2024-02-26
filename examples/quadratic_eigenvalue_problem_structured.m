%
% Test for the Backward error of a quadratic eigenvalue problem (QEP).
%

n = 3;
p = 2;

rng(5)
A0 = randn(n, n); A0 = -A0 * A0';
A1 = randn(n, n); A1 = A1 + A1';
A2 = randn(n, n); A2 = eye(n) + A2 * A2';

A0(1,[1 2])=zeros(1,2);
A0(3,[2 3])=zeros(1,2);
A0(2,[1 3])=zeros(1,2);

[V, L] = polyeig(A0, A1, A2);
f = @(x) [ 1, x, x^2 ];

% Select p eigenpairs to test
L = L(1:p); V = V(:, 1:p);

% Construct perturbations
epsilon = 1e-3;
dA0 = randn(n, n); dA0 = dA0 / norm(dA0) * epsilon;
dA1 = randn(n, n); dA1 = dA1 / norm(dA1) * epsilon;
dA2 = randn(n, n); dA2 = dA2 / norm(dA2) * epsilon;

dA0(1,[1 2])=zeros(1,2);
dA0(3,[2 3])=zeros(1,2);
dA0(2,[1 3])=zeros(1,2);

P1=zeros(n^2,3);
P1(3,1)=1;
P1(5,2)=1;
P1(7,3)=1;

P2=eye(n^2);
P3=eye(n^2);

P=blkdiag(P1,P2,P3);
%keyboard

% Setup for the problem
F = {A0 + dA0, A1 + dA1, A2 + dA2};

R = be_residual(F, f, V, diag(L));
fprintf('Norm of the residual for the perturbed problem: %e\n', norm(R, 'fro'));

D = be_unstructured(F, f, V, diag(L));
R = be_residual({F{1}+D{1}, F{2}+D{2}, F{3}+D{3}}, f, V, diag(L));
fprintf('Norm of the residual for the computed BE: %e\n', norm(R, 'fro'));

fprintf('Norm of the computed backward error: || D ||_F = %e\n', be_norm(D));

% Compute estimates for the backward error with the structured bound
    bnd = be_linear_structured_bound(F, f, V, diag(L),P);
    fprintf('STRUCT BOUND %e\n', bnd);