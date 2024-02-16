function Ft = be_perturb(F, epsilon)
%BE_PERTURB 

Ft = F;

w = randn(1, length(F));
w = w / norm(w) * epsilon;

for j = 1 : length(Ft)
    dF = randn(size(Ft{j}));
    dF = dF / norm(dF, 'fro');
    Ft{j} = F{j} + w(j) * dF;
end

end

