function [V, L] = be_newton(F, f, lstart, vstart)
%BE_NEWTON 

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
    N = inf;

    it = 0;

    while norm(N) > 1e-12
        it = it + 1;

        v = w(2:end); l = w(1);
        [fv, fpv] = f(l);

        RES = fv(1) * (F{1}*v);
        for jj = 2 : k
            RES = RES + fv(jj) * (F{jj} * v);
        end

        JFv = fpv(1) * (F{1} * v);
        for jj = 2 : k
            JFv = JFv + fpv(jj) * (F{jj} * v);
        end

        FF = fv(1) * F{1};
        for jj = 2 : k
            FF = FF + fv(jj) * F{jj};
        end

        RES = [ RES; e' * v - 1 ];
        J = [ JFv,  FF ; 0, e' ];
        
        N = J \ RES;
        % N(2:end) = conj(N(2:end));

        % ll = min(1, w(1) / (100 * N(1)));
        ll = 1;

        w = w - ll * N;

        % fprintf('Iteration %d, RES = %e, l = %e + %ei\n', it, norm(RES), real(l), imag(l));
        % keyboard;
    end

    V(:, j) = w(2:end) / norm(w(2:end));
    L(j,j)  = w(1);
end

end

