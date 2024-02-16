function nrm = be_norm(D)
%BE_NORM 

nrm = zeros(1, length(D));
for j = 1 : length(D)
    nrm(j) = norm(D{j}, 'fro');
end

nrm = norm(nrm);


end

