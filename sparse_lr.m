function A = sparse_lr(A, U, V)
%SPARSE_LR

[I, J] = find(A);
w = zeros(length(I), 1);

for s = 1 : length(I)
    w(s) = U(I(s), :) * V(J(s), :)';
end

A = sparse(I, J, w);

end

