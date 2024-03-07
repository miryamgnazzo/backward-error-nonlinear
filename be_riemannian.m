function [D] = be_riemannian(F, f, structures, V, L)
%BE_RIEMANNIAN
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
    n = size(F{1}{1}, 1);
else
    n = size(F{1}, 1);
end

% Build the manifold for the optimization problem
elements = struct;
for j = 1 : length(structures)
    MM = [];

    switch structures
        case 'low-rank'
            k = size(F{j}{1}, 2);
            MM = fixedrankembeddedfactory(n, n, k);
        case 'identity'
            % FIXME: We should implement the sparse factory here, and not
            % use the sparse one which does not force all elements equal
            % along the diagonal.
            MM = euclideansparsefactory(F{j});
        case 'sparse'
            MM = euclideansparsefactory(F{j});
    end

    elname = sprintf('F%d', j);
    elements.(elname) = MM;
end

M = productmanifold(elements);

% Select a starting point (we use the coefficients in F)
X0 = M.rand();
for j = 1 : length(structures)
    elname = sprintf('F%d', j);

    switch(structures)
        case 'low-rank'
            [QU, RU] = qr(F{j}{1}, 0); [QV, RV] = qr(F{j}{2}, 0);
            [U, S, V] = svd(RU * RV');
            X0.(elname) = struct('U', U * QU, 'S', S, 'V', V * QV);
        case 'sparse'
            XO.(elname) = F{j};
        case 'identity'
            XO.(elname) = F{j};
    end
end

X = X0;

% Select the regularization parameter, and repeat the optimization process
% while mu goes to zero.
mu = 1;
for j = 1 : 10
    X = be_riemannian_step(F, f, structures, mu, V, L, M, X);
    mu = mu / 4;
end

% Extract the perturbations from X
D = F;
for j = 1 : structures
    elname = sprintf('F%d', j);
    switch structures(j)
        case 'sparse'
            D{j} = X.(elname) - F{j};
        case 'identity'
            D{j} = X.(elname) - F{j};
        case 'low-rank'
            XL = X.(elname);
            D{j} = { [ XL.U * XL.S, -F{j}{1} ], [ XL.V, F{j}{2} ] };
    end
end
    
end

function X = be_riemannian_step(F, f, structures, mu, V, L, M, X0)

problem.M=M;
problem.cost = @(X) cost(F, f, structures, mu, V, L, X);

A=[A2 A1 A0];

% Euclidean gradient. Projection from ambient space to the tangent space of
% U(n) is handled automatically (see stiefelcomplexfactory documentation)
problem.egrad = @egrad;
%problem.grad = @(X) problem.M.egrad2rgrad(X, egrad(X));
% Euclidean Hessian. Projection is handled automatically.
%problem.ehess = @ehess;

%Y=M.randvec(X0);

%hessfd=getHessianFD(problem,X0,Y);
%keyboard

%Commenting the instruction for the gradient/Hessian we use the FD
%approximation in Manopt
%Uncommenting them we use the exact form provided in egrad or ehess

options.maxiter = maxiter;
options.maxtime = timemax;
options.tolgradnorm = 1e-12;
options.Delta_bar = 4.47214*1e-0;
options.Delta0 = options.Delta_bar/8;
options.debug=0;
options.rho_regularization = 1e3;

warning('on', 'manopt:getHessian:approx');

%START of the optimization (usually trustregions)
[X, xcost, info, options] = trustregions(problem,X0, options);
%[X, xcost, info, options] = rlbfgs(problem, V0, options);
%[X, xcost, info, options] = steepestdescent(problem, X0, options);

%keyboard


infotable = struct2table(info);
%e = sqrt(infotable.cost);
e = infotable.cost;
t = infotable.time;

%X will be the final output
disp('FINAL COST')
xcost


%Final pertubations from the output X
D2=X.F2;
DD1=X.F1;
D0=X.F0;

%From structure to matrix (for the low rank coefficient)
D1=(DD1.U)*(DD1.S)*(DD1.V)';

disp('FINAL Dist')
norm([D2 D1 D0],'fro')

end

function f = cost(F, f, structures, mu, V, L, X)
    RES = zeros(size(V));

    for j = 1 : size(V, 2)
        fv = f(L(j,j));

        for i = 1 : length(structures)
            elname = sprintf('F%d', i);
            switch structures(j)
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
    for j = 1 : structures
        elname = sprintf('F%d', j);
        switch structures(j)
            case 'low-rank'
                fr(j) = lr_norm(...
                    [ X.(elname).U * X.(elname).S, -F{j}{1} ], ...
                    [ X.(elname).V, F{j}{2} ]);
            case 'sparse'
                fr(j) = norm(X.(elname) - F{j}, 'fro');
            case 'identity'
                fr(j) = norm(X.(elname) - F{j}, 'fro');
        end
    end

    f = norm([ f, mu * fr ]);
end

function f = costold(X)
%Computes the cost function given an element X of the product manifold
    
   E2=X.F2;
   E1=X.F1;
   
   U=E1.U;
   V=E1.V;
   S=E1.S;
  
   E0=X.F0;
   
   EE1=U*S*V'; %From structure to matrix for the low rank term
   
   E=[E2 EE1 E0];
   W=[Vhat*L^2; Vhat*L; Vhat]; 
   
   f = mu*norm((A+E)*W,'fro')^2 + (norm(E,'fro')^2);     
%  f = norm((A+E)*W,'fro')^2 + mu*(norm(E,'fro')^2);     

end

function g = egrad(X)
%Computes the euclidian gradient given an element X of the product manifold

   E2 = X.F2;
   E1 = X.F1;
  
   U = E1.U;
   V = E1.V;
   S = E1.S;
   
   EE1=U*S*V'; %From structure to matrix for the low rank term

   E0 = X.F0;
      
   W2=Vhat*L^2; 
   W1=Vhat*L; 
   W0=Vhat;
   
   E = [E2 EE1 E0];
   W=[W2; W1; W0];
   
%   Matrix=F*W;
   Matrix=(E+A)*W;


   g.F0=2*Matrix*W0';
   g.F1=2*Matrix*W1'; 
   g.F2=2*Matrix*W2';
    
   g.F2=mu*g.F2 + 2*(E2);
   g.F1=mu*g.F1 + 2*(EE1);
   g.F0=mu*g.F0 + 2*(E0);

%    g.F2= g.F2 + 2*mu*(E2);
%    g.F1= g.F1 + 2*mu*(EE1);
%    g.F0= g.F0 + 2*mu*(E0);

   %In this way I usually compare the eucliedian gradient with the FD one
%   gradfun = approxgradientFD(problem);
%    temp1 = gradfun(X);
%    temp2 = problem.M.proj(X,g);
%        
%     nn=norm(temp1.F2 - temp2.F2,'f')/ norm(temp2.F2,'fro')
%     
%     keyboard
end

function nrm = lr_norm(U, V)
    [~, RU] = qr(U, 0);
    [~, RV] = qr(V, 0);
    nrm = norm(RU * RV', 'fro');
end