function [D] = be_unstructured_bound(bndtype, F, f, V, L)
%BE_UNSTRUCTURED Return the backward error D for the NEP. 
%
% F = { F1, ..., Fk } coefficient matrices
% f = { f1, ..., fk } functions
% V = [v1, ..., vp] approximate eigenvectors
% L = diag(l1, ...,lp) approximate eigenvalues
% bndtype may be 1,2,3
%
%
% Allowed BNDTYPE:
%  - 1: Bound based on the smallest singular value of F, the matrix with 
%       the evaluations of the functions at the eigenvalues. Only holds for
%       p <= k.
%  - 2: Bound based on the smallest singular value of F, and the condition 
%       number of V, the eigevector matrix. Only holds for p <= kn.
%  - 3: Bound based on smallest singular value of the Khatri-Rao product of
%       F and V.

p = size(V, 2);
k = length(F);
n = size(V, 1);

% For bound 2, we need the assumption that the columns of V are scaled to
% have unit norm; this is a good idea in any case, so we rescale them here
V = V ./ vecnorm(V);

R = be_residual(F, f, V, L);

% We handle the diagonal L case separately from the rest, because it can be
% considerably cheaper.
if isdiag(L)
    FF = zeros(p, k);
    for i = 1 : p
        fv = f(L(i,i));
        FF(i, :) = fv;
    end

    if bndtype == 1
        if p > k
            error('BNDTYPE 1 only holds for p <= k')
        end

        sp = min(svd(FF));
        D = norm(R, 'fro') / sp;
        return;
    end

    if bndtype == 2
        if p > k*n
            error('BNDTYPE 2 only holds for p <= nk')
        end
        sp = min(svd(FF));
        D = norm(R, 'fro') / sp * cond(V);
        return;
    end

    if bndtype == 3
        M = Khatri_Rao(FF, V');
        % We truncate under machine precision, to ignore zero singular
        % values that could be contaminated by roundoff errors.
        s = svd(M); s = s(s >= s(1) * eps);
        D = norm(R, 'fro') / min(s);
        return
    end

    error('Unsupported Bound Type requested')
else
    error('Not implemented')
end

end

function K=Khatri_Rao(F,V)
% K= Kathri- Rao transpose product between F and V

K=zeros(size(F,1),size(F,2)*size(V,2));

if (size(F,1)~=size(V,1))
    error('inconsistent dimension') 
end

for i=1:size(F,1)
    K(i,:)=kron(F(i,:),V(i,:)); 
end

end
