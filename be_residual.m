function R = be_residual(F, f, V, L)
%BE_RESIDUAL
% Residual R= F1 V f1(L) + .. + Fk V fk(L)
%
% F = { F1, ..., Fk } coefficient matrices
% f = { f1, ..., fk } functions
% V = [v1, ..., vp] approximate eigenvectors
% L = diag(l1,...,lp) approximate eigenvalues

if isempty(F)
    R = zeros(size(V));
    return;
end

if isdiag(L)
    R = zeros(size(V));
    for j = 1 : size(L, 2)
        fv = f(L(j,j));
        for i = 1 : length(F)
            R(:,j) = R(:,j) + fv(i) * (F{i} * V(:,j));
        end
    end
else
    error('Non-diagonal L are currently unsupported')
end


end

