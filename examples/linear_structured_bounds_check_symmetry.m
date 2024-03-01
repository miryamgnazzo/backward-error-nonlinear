% Symmetry case: comparison of the result with the structured error and the
% error bound for symmetric structure
% This script checks the structured bound that we have on randomly
% generated problems and compare this bound with the unstructured one
% (obtained with the Kathri Rao product).
%

ntests = 1000;
epsilon = 1e-3;
V = 0;
n = 64;
%p=2;

while cond(V) > 1e10
    F = cell(1, 5);

    F{1} = eye(n);
    F{2} = randn(n) / 10; F{2} = F{2} + F{2}';
    F{3} = randn(n); F{3} = F{3} + F{3}';
    F{4} = randn(n); F{4} = F{4} + F{4}';
    F{5} = randn(n); F{5} = F{5} + F{5}';
   
    [V, L] = be_newton(F, @f, -1 : 1);

end


fprintf('finito Newton \n');


 %Construction of the matrix P for the linear structure
    P1=symmetry_struct(ones(n));
    
    P=blkdiag(P1,P1,P1,P1,P1);

nrm = zeros(1, ntests);
be = zeros(1, ntests);
bnd = zeros(3, ntests);

for s = 1 : ntests
    % Ft = be_perturb(F, epsilon * exp(randn));
    Vt = V; Vt = Vt + randn(size(Vt)) * diag(epsilon * randn(1, size(Vt, 2)));
    Lt = L; Lt = Lt + diag(epsilon * randn(1, size(Lt, 1)));

    Ft = be_perturb(F, epsilon * exp(randn));

       %Eventualmente aggiungo symmetric structure
%     Ft{1} = Ft{1} + Ft{1}';
%     Ft{2} = Ft{2} + Ft{2}';
%     Ft{3} = Ft{3} + Ft{3}';
%     Ft{4} = Ft{4} + Ft{4}';
%     Ft{5} = Ft{5} + Ft{5}';

    R = be_residual(Ft, @f, V, L);
    nrm(s) = norm(R, 'fro');
    D = be_linear_structured(Ft, @f, V, L, P);

    be(s) = be_norm(D);

    %Unstructured bound 
    bnd(1, s) = be_unstructured_bound(3, Ft, @f, V, L);

    %Structured bound (linear structures)
    bnd(2, s) = be_linear_structured_bound(Ft, @f, V, L, P);

    %Bound for the symmetric case
    bnd(3, s) = be_symmetric_bound(Ft, @f, V, L);
    
end

% Sort them
[nrm, I] = sort(nrm);
be = be(I);
bnd = bnd(:, I);

figure;
loglog(nrm, be, 'r*'); hold on;
plot(nrm, bnd(1, :), 'k--');
plot(nrm, bnd(2, :), 'b--');
plot(nrm, bnd(3, :), 'g--');

writematrix([nrm', be', bnd'], './linear_structured_bounds_check_symm_bound.dat', 'Delimiter', '\t');

 function [fv, fvp] = f(x)
     fv = [ x^2, x, 1, expm(-x), expm(-2*x)];
    if nargout > 1
        fvp = [2*x, 1, 0, -expm(-x), -2*expm(-2*x)];
    end
end