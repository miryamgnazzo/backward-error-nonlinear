function [D] = be_unstructured(F, f, V, L)
%BE_UNSTRUCTURED
% Return the backward error D for the NEP. 
%
% F = { F1, ..., Fk }
% f = { f1, ..., fk }
% V = [v1, ..., vp] approximate eigenvectors
% L = diag(l1,...,lp) approximate eigenvalues

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

function K=Khatri_Rao(F,V)
%K = Kathi-Rao transpose product between F and V

K=zeros(size(F,1),size(F,2)*size(V,2));

if (size(F,1)~=size(V,1))
    error('inconsistent dimension') 
end

for i=1:size(F,1)
    K(i,:)=kron(F(i,:),V(i,:)); 
end

end