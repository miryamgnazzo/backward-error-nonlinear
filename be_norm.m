function nrm = be_norm(D)
%BE_NORM 

nrm = zeros(1, length(D));
for j = 1 : length(D)
    if isa(D{j}, 'lowrank')
        % low-rank
        nrm(j) = lr_norm(D{j}.U, D{j}.V, 'fro');
    else
        nrm(j) = norm(D{j}, 'fro');
    end
end

nrm = norm(nrm);


end

