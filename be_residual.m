function R = be_residual(F, f, V, L)
%BE_RESIDUAL

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
    error('Currently unsupported')
end


end

