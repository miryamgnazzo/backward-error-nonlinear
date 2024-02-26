function [D] = be_linear_structured(F, f, V, L, P)
%BE_LINEAR_STRUCTURED Return the backward error D for the NEP. 
%
% F = { F1, ..., Fk }
% f = { f1, ..., fk }
% P matrix is a blkdiag matrix containing the information on the structure
% of F


p = size(V, 2);
k = length(F);
n = size(V, 1);

R = be_residual(F, f, V, L);
r = reshape(R,n*p,1);

% We handle the diagonal L case separately from the rest, because it can be
% considerably cheaper.
if isdiag(L)
    FF = zeros(p, k);
    FF_t = zeros(p,p*k);

    for i = 1 : p
        fv = f(L(i,i));
        FF(i, :) = fv;
    end


   %QUESTO SI Può fare meglio sicuramente
    for j=1:k
      FF_t(:,1+(j-1)*p:j*p)=diag(FF(:,j));
    end

    %bound for linear structures
    % We create the matrix M
    M=FF_t*kron(eye(k),V');
    M=kron(M,eye(n))*P;
    

   % M = Khatri_Rao(FF, V');
   % DD = -R * pinv(M');
    
    delta = -pinv(M)*r;
   % DD = P*delta;
    
    DD = reshape(P*delta,n,n*k);

    D = cell(1, k);
    for j = 1 : k
        D{j} = DD(:, (j-1)*n+1:j*n);
    end
else
    error('Not implemented')
end

end