function [D] = be_symmetric(F, f, V, L)
%BE_SYMMETR Return the backward error D for the NEP, symmetric case
%
% F = { F1, ..., Fk } coefficient matrices
% f = { f1, ..., fk } functions
% V = [v1, ..., vp] approximate eigenvectors
% L = diag(l1, ..., lp) approximate eigenvalues

p = size(V, 2);
k = length(F);
n = size(V, 1);

Res = be_residual(F, f, V, L);

% We handle the diagonal L case separately from the rest, because it can be
% considerably cheaper.
if isdiag(L)
   
    FF = zeros(p, k);

    for i = 1 : p
        fv = f(L(i,i));
        FF(i, :) = fv;
    end
 
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

    % Solve for Blocks 21
    Block_21 = Res2 * pinv(T);

    F21t = cell(1, k);

    for m = 1:k
        F21t{m} = Block_21(:, (m-1)*p + 1 : m*p);
    end

    % Solve for the symmetric blocks F11
    P = 1 : p^2; P = reshape(P, p, p); P = P'; P = P(:);
    CM = eye(p^2); CM = CM(:, P);
    
    Blk = zeros(k*p^2);

    for h = 1 : k
        Blk((h-1)*p^2 + 1 : h*p^2 , (h-1)*p^2 + 1 : h*p^2 ) = CM - eye(p^2);
    end

    Sys_11 = [kron( T' , eye(p)) ; Blk];


    Block_11 = pinv(Sys_11) * [ Res1(:) ; zeros(k*p^2, 1) ];
    
    F11t = cell(1, k);

    for q = 1 : k
        F11t{q} = reshape(Block_11((q-1)*p^2 + 1 : q*p^2), p, p);
    end

    %Final solution
    D = cell(1, k);

    for j = 1 : k
        D{j} = [ F11t{j}, F21t{j}' ; F21t{j} , zeros(n-p, n-p) ];
        D{j} = Q * D{j} * Q';
    end

else
    error('Not implemented')
end

end