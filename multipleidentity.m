function M = multipleidentity(n)
    % n = size of the matrices involved
    
    
    dimensions_vec = n*ones(2,1);
    assert(length(dimensions_vec) == 2, 'A should be a matrix (or a vector).');
    S = speye(n);
    
    M.size = @() dimensions_vec;
    
    M.name = @() sprintf('Euclidean space R^(%dx%d) with multiple of identity', ...
                                            dimensions_vec(1), dimensions_vec(2));
    
    M.dim = @() 1;

    M.inner = @(x, d1, d2) dimensions_vec(1)*d1(1,1)'*d2(1,1); %d1(:).'*d2(:); 
    
    M.norm = @(x, d) sqrt(dimensions_vec(1))*abs(d(1,1)); %norm(d, 'fro');
    
    M.dist = @(x, y) sqrt(dimensions_vec(1))*abs(x(1,1)-y(1,1));

    M.typicaldist = @() sqrt(prod(dimensions_vec));
    
    M.proj = @(x, d) S.*d; 


end