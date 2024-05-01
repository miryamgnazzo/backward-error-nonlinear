function [D] = be_linear_structured_bound(F, f, V, L,P)
%BE_LINEAR_STRUCTURED_BOUND Return the backward error D for the NEP.
%
% F = { F1, ..., Fk }
% f = { f1, ..., fk }
% P = block diagonal matrix with the information on the structure
%     separetely on the coefficients P^1,...,P^k (these are not square) 
%      already with orthonormal columns in P^i
%
%
% D = Bound based on smallest singular value of the matrix
% ((FF \odot^T V^T)\otimes I_n) P

p = size(V, 2);
k = length(F);
n = size(V, 1);

% For bound 2 (unstructured), we need the assumption that the columns of V are scaled to
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

   % We create the matrix M
   M = Khatri_Rao(FF, V');
   M=kron(M,eye(n))*P;
    
    % We truncate under machine precision, to ignore zero singular
    % values that could be contaminated by roundoff errors.
    s = svd(M); s = s(s >= s(1) * eps);
    D = norm(R, 'fro') / min(s);

else
    error('Not implemented')
end

end

function K=Khatri_Rao(F,V)
% Khatri Rao transpose between F and V
    K=zeros(size(F,1),size(F,2)*size(V,2));
    
    if (size(F,1)~=size(V,1))
        error('inconsistent dimension') 
    end
    
    for i=1:size(F,1)
        K(i,:)=kron(F(i,:),V(i,:)); 
    end
end