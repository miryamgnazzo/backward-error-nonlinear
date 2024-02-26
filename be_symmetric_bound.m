function [D] = be_symmetric_bound(F, f, V, L)
%BE_SYMMETRIC_BOUND Return the backward error D for the NEP.
%
% F = { F1, ..., Fk }
% f = { f1, ..., fk }
%
%
% D=Bound for the symmetric case

p = size(V, 2);
k = length(F);
n = size(V, 1);

% For bound 2, we need the assumption that the columns of V are scaled to
% have unit norm; this is a good idea in any case, so we rescale them here
V = V ./ vecnorm(V);

Res = be_residual(F, f, V, L);

% We handle the diagonal L case separately from the rest, because it can be
% considerably cheaper.
if isdiag(L)
    FF = zeros(p, k);

    for i = 1 : p
        fv = f(L(i,i));
        FF(i, :) = fv;
    end 


    %Construction of the matrix \tilde T (here =T)
    [Q, R] = qr(V);
    Qt  = Q(:, 1:p);
    Rt = R(1:p, :);
    Qto = Q(:, p+1:end);
    
    Res1 = - Qt' * Res;
    Res2 = - Qto' * Res;

    T = zeros(k*p,p);

    for  l = 1 : k 
        T((l-1)*p + 1: l*p,:) =  Rt * diag(FF(:, l));
    end

    % Construction matrix Ms (here is Sys_11)
    P = 1 : p^2; P = reshape(P, p, p); P = P'; P = P(:);
    CM = eye(p^2); CM = CM(:, P);
    
    Blk = zeros(k*p^2);

    for h = 1 : k
        Blk((h-1)*p^2 + 1 : h*p^2 , (h-1)*p^2 + 1 : h*p^2 ) = CM - eye(p^2);
    end

    Sys_11 = [kron( T' , eye(p)) ; Blk];

   bound = norm(pinv(Sys_11),'fro')^2 + 2*norm(T,'fro')^2;
   D= norm(Res,'fro')*sqrt(bound);

else
    error('Not implemented')
end

end