%
% This script checks the unstructured bounds that we have on "beam"
% problem from NEP-PACK package.
%
% Return a comparison among the possible upper bounds for 1000 tests on the
% beam problem, with perturberd coefficients Ft.

n = 1000;
ntests = 1000;
epsilon = 1e-3;
V = 0;
p = 10;

A0 = spdiags(ones(n, 1) * [1 -2 1], -1:1, n, n);
A0(end, end-1:end) = [ -n, n ];
A1 = sparse(n,n,1);
F = { speye(n), A0, A1 };
[V, L] = be_newton(F, @f, -linspace(0.1, 1, p));

nrm = zeros(1, ntests);
be = zeros(1, ntests);
bnd = zeros(3, ntests);

for s = 1 : ntests

    Ft = be_perturb(F, epsilon * exp(randn));
    R = be_residual(Ft, @f, V, L);
    nrm(s) = norm(R, 'fro');
    D = be_unstructured(Ft, @f, V, L);

    be(s) = be_norm(D);

    if (p <= 3)

        for j = 1 : 3
           bnd(j, s) = be_unstructured_bound(j, Ft, @f, V, L);
        end

    else
    %BNDTYPE = 1 not support for p > 3 
        for j = 2 : 3
            bnd(j, s) = be_unstructured_bound(j, Ft, @f, V, L);
        end

    end
end

% Sort them
[nrm, I] = sort(nrm);
be = be(I);
bnd = bnd(:, I);

figure;
if (p <= 3)
    loglog(nrm, be, 'r*'); hold on;
    plot(nrm, bnd(1, :), 'k--');
    plot(nrm, bnd(2, :), 'b--');
    plot(nrm, bnd(3, :), 'g--');
else
% BNDTYPE = 1 not supported for p > 3    
    loglog(nrm, be, 'r*'); hold on;
    plot(nrm, bnd(2, :), 'b--');
    plot(nrm, bnd(3, :), 'g--');
end

%writematrix([nrm', be', bnd'], './unstructured_bounds_check_beam.dat', 'Delimiter', '\t');

function [fv, fp] = f(lambda)
    fv = [ -lambda, 1.0,  exp(-lambda) ];
    fp = [ -1.0,    0.0, -exp(-lambda) ];
end