function nrm = be_norm(D)
%BE_NORM 
% Norm of the cell array D
%
% D = {D1, ..., Dk} cell array

nrm = zeros(1, length(D));
for j = 1 : length(D)
    if isa(D{j}, 'lowrank')
       % low-rank
       nrm(j) = lr_norm(D{j}.U, D{j}.V);
    else
        nrm(j) = norm(D{j}, 'fro');
    end
end

nrm = norm(nrm);

end

function nrm = lr_norm(U, V)
%Computation of the norm for low rank matrices UV'
    [~, RU] = qr(U, 0);
    [~, RV] = qr(V, 0);
    nrm = norm(RU * RV', 'fro');
end