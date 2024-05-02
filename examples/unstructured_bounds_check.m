%
% This script checks the unstructured bounds that we have on randomly
% generated problems.
%

ntests = 1000;
epsilon = 1e-3;
V = 0;
n = 128;

while cond(V) > 1e10
    F = cell(1, 5);
    F{1} = eye(n);
    F{2} = randn(n) / 10; F{2} = F{2} * F{2}';
    F{3} = randn(n); F{3} = -F{3} * F{3}';
    F{4} = randn(n); F{4} = F{4} * F{4}';
    F{5} = randn(n); F{5} = F{5} * F{5}';
    
    % approximation of 3 eigenpairs
    [V, L] = be_newton(F, @f, -1 : 1);
end

nrm = zeros(1, ntests);
be = zeros(1, ntests);
bnd = zeros(3, ntests);

for s = 1 : ntests
    Vt = V; Vt = Vt + randn(size(Vt)) * diag(epsilon * randn(1, size(Vt, 2)));
    Lt = L; Lt = Lt + diag(epsilon * randn(1, size(Lt, 1)));

    Ft = be_perturb(F, epsilon * exp(randn));
    R = be_residual(Ft, @f, V, L);
    nrm(s) = norm(R, 'fro');
    D = be_unstructured(Ft, @f, V, L);

    be(s) = be_norm(D);
    for j = 1 : 3
        bnd(j, s) = be_unstructured_bound(j, Ft, @f, V, L);
    end
end

% Sort them
[nrm, I] = sort(nrm);
be = be(I);
bnd = bnd(:, I);

figure;
loglog(nrm, be, 'r*'); hold on;
plot(nrm, bnd(1, :), 'k--');
plot(nrm, bnd(2, :), 'b--');
plot(nrm, bnd(3, :), 'm--');

function [fv, fvp] = f(x)
    fv = [ x^2, x, 1, expm(-x), expm(-2*x)];
    if nargout > 1
        fvp = [2*x, 1, 0, -expm(-x), -2*expm(-2*x)];
    end
end
