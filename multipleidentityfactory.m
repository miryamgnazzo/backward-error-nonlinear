function M = multipleidentityfactory(n)
    % n = size of the matrices involved
      
    M.size = @() [n n];
    
    M.name = @() sprintf('Euclidean space R^(%dx%d) with multiple of identity', ...
                                            n, n);
    
    M.dim = @() 1;

    M.inner = @(x, d1, d2) n*d1'*d2;
    
    M.norm = @(x, d) sqrt(n)*abs(d); %norm(d, 'fro');
    
    M.dist = @(x, y) sqrt(n) * abs(x - y);

    M.typicaldist = @() sqrt(n);
    
    M.proj = @proj;

    function l = proj(x, d)
        if isscalar(d)
            l = d;
        else
            l = trace(d) / n;
        end
    end

    M.egrad2rgrad = @(x, g) trace(g) / n;
    
    M.ehess2rhess = @(x, eg, eh, d) trace(eh) / n;
    
    M.tangent = M.proj;
    
    M.exp = @exp;
    function y = exp(x, d, t)
        if nargin == 3
            y = x + t * d;
        else
            y = x + d;
        end
    end
    
    M.retr = M.exp;
    
    M.log = @(x, y) y-x;

    M.hash = @(x) sprintf('%1.16f', x);
    
    M.rand = @() randn;
    
    M.randvec = @(x) sign(rand - .5) / sqrt(n);
    
    M.lincomb = @lincomb;

    function v = lincomb(x, a1, u1, a2, u2)
        if ~exist('a2', 'var')
            a2 = 0;
            u2 = 0;
        end
        v = a1 * u1 + a2 * u2;
    end
    
    M.zerovec = @(x) 0;
    
    M.transp = @(x1, x2, d) d;

    M.isotransp = M.transp; % the transport is isometric
    
    M.pairmean = @(x1, x2) .5*(x1+x2);
    
    M.vec = @(x, u_mat) u_mat;
    M.mat = @(x, u_vec) u_vec;
    M.vecmatareisometries = @() true;

    M.lie_identity = @() 1;
end