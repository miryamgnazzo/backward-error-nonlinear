%
% This script checks the unstructured bounds that we have on Hadeler
% function in nlevp package.
%


% Hadeler problem from the NLEVP collection
%

ntests = 1000;
epsilon = 1e-3;
V = 0;

n = 8;
[F, f]=nlevp('hadeler', n);
%keyboard

% The first two eigenvalues are good initial estimate for n = 8
[V, L] = be_newton(F, f, [0.217, 0.885, 4]);

nrm = zeros(1, ntests);
be = zeros(1, ntests);
bnd = zeros(3, ntests);

for s = 1 : ntests
%    Vt = V; Vt = Vt + randn(size(Vt)) * diag(epsilon * randn(1, size(Vt, 2)));
%    Lt = L; Lt = Lt + diag(epsilon * randn(1, size(Lt, 1)));

    Ft = be_perturb(F, epsilon * exp(randn));
    R = be_residual(Ft, f, V, L);
    nrm(s) = norm(R, 'fro');
    D = be_unstructured(Ft, f, V, L);

    be(s) = be_norm(D);
    for j = 1 : 3
        bnd(j, s) = be_unstructured_bound(j, Ft, f, V, L);
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

writematrix([nrm', be', bnd'], './unstructured_bounds_check_hadeler.dat', 'Delimiter', '\t');

% function [fv, fvp] = f(x)
%     fv = [ x^2, x, 1, expm(-x), expm(-2*x)];
%     if nargout > 1
%         fvp = [2*x, 1, 0, -expm(-x), -2*expm(-2*x)];
%     end
% end