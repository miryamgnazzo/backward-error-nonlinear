function [V, L] = be_newton(F, f, lstart)
%BE_NEWTON 

p = length(lstart);
n = size(F{1}, 1);
V = zeros(n, p);
L = zeros(p);
k = length(F);

for j = 1 : p
    v = eye(n, 1);
    l = lstart(j);
    
    e = ones(n, 1);
    w = [ l ; v ];    
    N = inf;

    while norm(N) > 1e-12
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
        w = w - N;
    end

    V(:, j) = w(2:end);
    L(j,j)  = w(1);
end

end

