function [D] = be_riemannian_hessian(F, f, structures, V, L)
%BE_RIEMANNIAN
% WITH the addition of the Riemannian Hessian for the cost function
%
% STRUCTURES should be a cell-array: { S1, ..., Sk }, where each Sj is one
% of the following:
%
%  - 'sparse'
%  - 'low-rank'
%  - 'identity': multiple of the identity matrix
%
% The coefficients should be specificed in structured form:
%
% - Sparse matrices should be stored in sparse MATLAB format.
%
% - low-rank coefficients should be in factored form, by passing a
%   cell-array with two matrices { U, V } such that F{j} = U*V'
%
% - multiple of the identities should be passed in sparse MATLAB format.
%
% Example: if A0, A1, A2 have structures 'identity', 'low-rank', and
% 'sparse', one should use
%
%  F = { l*speye(n), { U, V }, sparse(...) };

% Determine the size n of the problem (for low-rank problems, we look at
% the number of rows in U).
if strcmp(structures{1}, 'low-rank')
    n = size(F{1}.U, 1);
else
    n = size(F{1}, 1);
end

% Build the manifold for the optimization problem
elements = struct;
for j = 1 : length(structures)
    MM = [];

    switch structures{j}
        case 'low-rank'
            k = size(F{j}.U, 2);
            MM = fixedrankembeddedfactory(n, n, k);
        case 'identity'
            % FIXME: We should implement the identity factory here, and not
            % use the sparse one which does not force all elements equal
            % along the diagonal.
            MM = multipleidentityfactory(size(F{j}, 1));
        case 'sparse'
            MM = euclideansparsefactory(F{j});
    end

    elname = sprintf('F%d', j);
    elements.(elname) = MM;
end

M = productmanifold(elements);

% Select a starting point (we use the coefficients in F)
X = struct;
for j = 1 : length(structures)
    elname = sprintf('F%d', j);

    switch structures{j}
        case 'low-rank'
            [QU, RU] = qr(F{j}.U, 0); [QV, RV] = qr(F{j}.V, 0);
            [UU, S, VV] = svd(RU * RV');
            X.(elname) = struct('U', QU * UU, 'S', S, 'V', QV * VV);
        case 'sparse'
            X.(elname) = F{j};
        case 'identity'
            X.(elname) = elements.(elname).proj([], F{j});
    end
end

% Select the regularization parameter, and repeat the optimization process
% while mu goes to zero.
mu = 1;
while sqrt(mu) > eps
    fprintf('Running Riemannian optimization with mu = %e\n', mu);
    X = be_riemannian_step(F, f, structures, mu, V, L, M, elements, X);    
    mu = mu / 64;
end
% X = be_riemannian_step(F, f, structures, eps^2, V, L, M, elements, X);

% Extract the perturbations from X
D = F;
for j = 1 : length(structures)
    elname = sprintf('F%d', j);
    switch structures{j}
        case 'sparse'
            D{j} = X.(elname) - F{j};
        case 'identity'
            D{j} = X.(elname) * speye(n) - F{j};
        case 'low-rank'
            XL = X.(elname);
            D{j} = lowrank([ XL.U * XL.S, -F{j}.U ], [ XL.V, F{j}.V ]);
    end
end
    
end

function X = be_riemannian_step(F, f, structures, mu, V, L, M, elements, X0)

problem.M = M;
problem.cost = @(X) cost(F, f, structures, mu, V, L, X);
problem.grad = @(X) grad(F, f, structures, mu, V, L, elements, X);
problem.hess = @(X, dX) hess(F, f, structures, mu, V, L, elements, X, dX);

% keyboard;

options.maxiter = 2 + (mu == 1 || mu == 0) * 18;
options.maxtime = inf;
options.tolgradnorm = sqrt(mu);
options.Delta_bar = 4.47214*1e-0;
options.Delta0 = options.Delta_bar/8;
options.debug = 0;
options.rho_regularization = 1e3;
options.maxinner = 2000;

warning('on', 'manopt:getHessian:approx');

%START of the optimization (usually trustregions)
[X, xcost, info, options] = trustregions(problem, X0, options);
%[X, xcost, info, options] = rlbfgs(problem, X0, options);
%[X, xcost, info, options] = steepestdescent(problem, X0, options);

infotable = struct2table(info);
e = infotable.cost;
t = infotable.time;

end

function f = cost(F, f, structures, mu, V, L, X)
    RES = zeros(size(V));

    for j = 1 : size(V, 2)
        fv = f(L(j,j));

        for i = 1 : length(structures)
            elname = sprintf('F%d', i);
            switch structures{i}
                case 'sparse'
                    RES(:,j) = RES(:,j) + fv(i) * (X.(elname) * V(:,j));
                case 'identity'
                    RES(:,j) = RES(:,j) + fv(i) * (X.(elname) * V(:,j));
                case 'low-rank'
                    RES(:,j) = RES(:,j) + ...
                        fv(i) * X.(elname).U * X.(elname).S ...
                            * (X.(elname).V' * V(:,j));
                
            end
        end
    end

    f = norm(RES, 'fro');

    % Regularization term mu * dist(X, X0)
    fr = zeros(1, length(structures));
    for j = 1 : length(structures)
        elname = sprintf('F%d', j);
        switch structures{j}
            case 'low-rank'
                fr(j) = lr_norm(...
                    [ X.(elname).U * X.(elname).S, -F{j}.U ], ...
                    [ X.(elname).V, F{j}.V ]);
            case 'sparse'
                fr(j) = norm(X.(elname) - F{j}, 'fro');
            case 'identity'
                fr(j) = norm(X.(elname) - diag(F{j}), 2);
        end
    end

    % f = f^2 + mu * sum(fr.^2);
    f = norm([ f, sqrt(mu) * fr ]).^2;
end

function g = grad(F, f, structures, mu, V, L, elements, X)
    % Working assumption: L is diagonal, f_j(L) = diag(f_j(L_{ii}))
    % We may generalize this, but we need a "matrix-function" handle for
    % f_j, which is not available for NLEVP problems. 

    k = length(structures);
    n = size(V, 1);

    % Compute the block vector W = [ V f_1(L) ; ... ; V f_k(L) ]
    W = kron(ones(k, 1), V);
    for j = 1 : size(W, 2)
        l = L(j,j);
        fv = f(l);
        for i = 1 : length(fv)
            % This is V(:,j) * f_i(L(j,j))
            W((i-1)*n+1 : i*n, j) = W((i-1)*n+1 : i*n, j) * fv(i);
        end
    end

    % Compute Ft * W = F1 * W1 + ... + Fk * Wk;
    FtW = zeros(n, size(W, 2));
    for i = 1 : k
        elname = sprintf('F%d', i);
        Wi = W((i-1)*n+1:i*n, :);
        switch structures{i}
            case 'low-rank'
                FtW = FtW + X.(elname).U * X.(elname).S * (X.(elname).V' * Wi);
            case 'sparse'
                FtW = FtW + X.(elname) * Wi;
            case 'identity'
                FtW = FtW + X.(elname) * Wi;
        end
    end

    % Assemble all components of the gradient
    g = struct;
    for i = 1 : k
        % grad_i is 2 * (FtW * Wi' + mu * (X.Fi - Fi))
        elname = sprintf('F%d', i);
        Wi = W((i-1)*n+1:i*n, :);

        switch structures{i}
            case 'low-rank'
                egrad_i = struct;
                egrad_i.U = [ 2 * FtW, ...
                    2*mu*X.(elname).U, ...
                    -2*mu*F{i}.U ];
                egrad_i.S = blkdiag(eye(size(W, 2)), X.(elname).S, eye(size(F{i}.U, 2)));
                egrad_i.V = [ Wi, X.(elname).V, F{i}.V ];

                grad_i = elements.(elname).egrad2rgrad(X.(elname), egrad_i);
            case 'sparse'
                grad_i = 2 * sparse_lr(F{i}, FtW, Wi) + ...
                    2 * mu * (X.(elname) - F{i});
            case 'identity'
                grad_i = 2 * trace(Wi' * FtW) / n + ...
                    2 * mu * (X.(elname) - trace(F{i}) / n);
        end

        g.(elname) = grad_i;
    end
end

function h = hess(F, f, structures, mu, V, L, elements, X, dX)

    % keyboard

    % Working assumption: L is diagonal, f_j(L) = diag(f_j(L_{ii}))
    % We may generalize this, but we need a "matrix-function" handle for
    % f_j, which is not available for NLEVP problems. 
    
    k = length(structures);
    n = size(V, 1);

    % Compute the block vector W = [ V f_1(L) ; ... ; V f_k(L) ]
    W = kron(ones(k, 1), V);
    for j = 1 : size(W, 2)
        l = L(j,j);
        fv = f(l);
        for i = 1 : length(fv)
            % This is V(:,j) * f_i(L(j,j))
            W((i-1)*n+1 : i*n, j) = W((i-1)*n+1 : i*n, j) * fv(i);
        end
    end
    
    dXtemp = dX;

    % Compute Ft * W = F1 * W1 + ... + Fk * Wk;
    % Compute Et * W = E1 * W1 + ... + Ek * Wk; with dX.(elname)=Ei
    FtW = zeros(n, size(W, 2));
    EtW = zeros(n, size(W, 2));
    for i = 1 : k
        elname = sprintf('F%d', i);
        Wi = W((i-1)*n+1:i*n, :);
        switch structures{i}
            case 'low-rank'
               %the computation of the Euclidean hessian requires to
               %consider vectors in the ambient space, using M.tangent2ambient(x, u)
               dXtemp.(elname) = elements.(elname).tangent2ambient(X.(elname), dX.(elname));
               EtW = EtW + dXtemp.(elname).U * dXtemp.(elname).S * (dXtemp.(elname).V' * Wi);
               FtW = FtW + X.(elname).U * X.(elname).S * (X.(elname).V' * Wi);
            case 'sparse'
                EtW = EtW + dX.(elname) * Wi;
                FtW = FtW + X.(elname) * Wi;
            case 'identity'
                EtW = EtW + dX.(elname) * Wi;
                FtW = FtW + X.(elname) * Wi;
        end
    end

    % Assemble all components of the hessian
    h = struct;
    for i = 1 : k
        % hess_i is 2 *( EtW * W_i' + mu* dX_i) 
        % grad_i is 2 * (FtW * Wi' + mu * (X.Fi - Fi))
        elname = sprintf('F%d', i);
        Wi = W((i-1)*n+1:i*n, :);

        switch structures{i}
            case 'low-rank'    
                %Euclidean gradient needed for ehess2rhess
                egrad_i = struct;
                egrad_i.U = [ 2 * FtW, ...
                    2*mu*X.(elname).U, ...
                    -2*mu*F{i}.U ];
                egrad_i.S = blkdiag(eye(size(W, 2)), X.(elname).S, eye(size(F{i}.U, 2)));
                egrad_i.V = [ Wi, X.(elname).V, F{i}.V ];

               %Euclidean Hessian needed for ehess2rhess
               ehess_i = struct;             
               ehess_i.U = [ 2 * EtW, ...
                             2* mu *dXtemp.(elname).U ];
               ehess_i.S = blkdiag(eye(size(W, 2)), dXtemp.(elname).S);
               ehess_i.V = [ Wi, dXtemp.(elname).V];
                           
               hess_i = elements.(elname).ehess2rhess(X.(elname), egrad_i, ehess_i, dX.(elname));
           case 'sparse'
                hess_i = 2 * sparse_lr(F{i}, EtW, Wi) + ...
                    2 * mu * dX.(elname);
           case 'identity'
                hess_i = 2 * trace(Wi' * EtW) / n + ...
                    2 * mu * dX.(elname);
                
        end
        
        h.(elname) = hess_i;
    end
end


function nrm = lr_norm(U, V)
    [~, RU] = qr(U, 0);
    [~, RV] = qr(V, 0);
    nrm = norm(RU * RV', 'fro');
end