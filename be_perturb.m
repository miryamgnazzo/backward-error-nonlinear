function Ft = be_perturb(F, epsilon)
%BE_PERTURB 
% Perturbation of the coefficient matrices of size eps
%
% F = { F1, ..., Fk }
% eps = perturbation size

Ft = F;

w = randn(1, length(F));
w = w / norm(w) * epsilon;

for j = 1 : length(Ft)
    dF = randn(size(Ft{j}));
    dF = dF / norm(dF, 'fro');
    Ft{j} = F{j} + w(j) * dF;
end

end

