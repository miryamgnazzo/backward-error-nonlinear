function [D] = be_unstructured(F, f, V, L)
%BE_UNSTRUCTURED Return the backward error D for the NEP. 
%
% F = { F1, ..., Fk }
% f = { f1, ..., fk }

p = size(V, 2);
k = length(F);
n = size(V, 1);

R = be_residual(F, f, V, L);

% We handle the diagonal L case separately from the rest, because it can be
% considerably cheaper.
if isdiag(L)
    FF = zeros(p, k);
    for i = 1 : p
        fv = f(L(i,i));
        FF(i, :) = fv;
    end

    M = Khatri_Rao(FF, V');
    DD = -R * pinv(M');

    D = cell(1, k);
    for j = 1 : k
        D{j} = DD(:, (j-1)*n+1:j*n);
    end
else
    error('Not implemented')
end

end

