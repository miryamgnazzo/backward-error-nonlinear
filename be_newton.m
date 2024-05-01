function [V, L] = be_newton(F, f, lstart, vstart)
%BE_NEWTON 
% Newton method for the numerical approximation of p eigenpairs of
% nonlinear eigenvalue problems
%
%
% F = { F1, ..., Fk } coefficient matrices
% f = { f1, ..., fk } functions
% lstart = [l1, ..., lp] initial approximate eigenvalues
% vstart = [v1, ..., vp] initial approximate eigenvectors (optional)

p = length(lstart);
n = size(F{1}, 1);
V = zeros(n, p);
L = zeros(p, p);
k = length(F);
e = randn(n, 1);

for j = 1 : p

    if exist('vstart', 'var')
        v = vstart(:, p);
    else
        v = randn(n, 1);
    end
    v = v / (e' * v);

    l = lstart(j);
    
    w = [ l ; v ];    
    RES = inf;

    it = 0;

    while norm(RES) > 1e-12
        it = it + 1;

        v = w(1:end-1); l = w(end);
        [fv, fpv] = f(l);

        RES = fv(1) * (F{1}*v);
        for jj = 2 : k
            RES = RES + fv(jj) * (F{jj} * v);
        end

        % fprintf('Newton, step %d, RES norm = %e\n', it, norm(RES));

        JFv = fpv(1) * (F{1} * v);
        for jj = 2 : k
            JFv = JFv + fpv(jj) * (F{jj} * v);
        end

        FF = fv(1) * F{1};
        for jj = 2 : k
            FF = FF + F{jj} * fv(jj);
        end

        RES = [ RES; e' * v - 1 ];
        % J = [ JFv,  FF ; 0, e' ];
        
        if issparse(FF)
            % Solution through Schur complementation to exploit sparsity
            [ii, jj, vv] = find(FF);
            ii = [ ii ; (1 : n)' ; (n+1)*ones(n,1) ]; 
            jj = [ jj ; (n+1)*ones(n,1) ; (1:n)' ];
            vv = [ vv ; JFv ; e ];
            J = sparse(ii, jj, vv, n+1, n+1);

            N = J \ RES;
        else
            J = [ FF , JFv ; e', 0 ];
            N = J \ RES;
        end
        ll = 1;

        w = w - ll * N;

        % fprintf('Iteration %d, RES = %e, l = %e + %ei\n', it, norm(RES), real(l), imag(l));
    end

    V(:, j) = w(1:end-1) / norm(w(1:end-1));
    L(j,j)  = w(end);
end

end

