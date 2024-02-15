function [D0, D1, D2] = sparse_QEP(A0, A1, A2, V, L)

maxiter = 10000;
timemax = 100;

elements = struct();
elements.F2 = euclideansparsefactory(A2); %coefficient of degree 2
elements.F1=euclideansparsefactory(A1); %coefficient of degree 1
elements.F0 = euclideansparsefactory(A0); %coefficient of degree 0
M = productmanifold(elements);

problem.M=M;
problem.cost = @cost;

problem.grad = @grad;

% Precompute the Hessian, which is constant
W2 = V*L^2;
W1 = V*L;
W0 = V;
problem.hess = @(X, U) struct(...
    'F0', 2 * (sparse_lr(A0, U.F0 * W0, W0) + sparse_lr(A0, U.F1 * W1, W0) + sparse_lr(A0, U.F2 * W2, W0)), ...
    'F1', 2 * (sparse_lr(A1, U.F0 * W0, W1) + sparse_lr(A1, U.F1 * W1, W1) + sparse_lr(A1, U.F2 * W2, W1)), ...
    'F2', 2 * (sparse_lr(A2, U.F0 * W0, W2) + sparse_lr(A2, U.F1 * W1, W2) + sparse_lr(A2, U.F2 * W2, W2)));

% problem.egrad = @grad
% problem.ehess = @(X, U) struct(...
%     'F0', 2 * (U.F0 * W0 * W0' + U.F1 * W1 * W0' + U.F2 * W2 * W0'), ...
%     'F1', 2 * (U.F0 * W0 * W1' + U.F1 * W1 * W1' + U.F2 * W2 * W1'), ...
%     'F2', 2 * (U.F0 * W0 * W2' + U.F1 * W1 * W2' + U.F2 * W2 * W2'));

options.maxiter = maxiter;
options.maxtime = timemax;
options.tolgradnorm = 1e-6;
%options.Delta_bar = 4.47214*1e-0;
%options.Delta0 = options.Delta_bar/8;
options.debug=0;
options.rho_regularization = 1e3;

warning('on', 'manopt:getHessian:approx');

% Initial point
X0=M.rand();

X0.F2=A2;
X0.F1=A1;
X0.F0=A0;

% [U0, U1, U2] = unstruct_error(full(A0), full(A1), full(A2), V, L);
% X0.F2 = X0.F2 + (X0.F2 ~= 0) .* U2;
% X0.F1 = X0.F1 + (X0.F1 ~= 0) .* U1;
% X0.F0 = X0.F0 + (X0.F0 ~= 0) .* U0;

% checkgradient(problem);
% keyboard

% checkhessian(problem);
% keyboard;


%START of the optimization (usually trustregions)
[X, xcost, info, options] = trustregions(problem, X0, options);
%[X, xcost, info, options] = rlbfgs(problem, V0, options);
%[X, xcost, info, options] = steepestdescent(problem, X0, options);

%keyboard
infotable = struct2table(info);
e = sqrt(infotable.cost);
t = infotable.time;

%X will be the final output
disp('FINAL COST')
xcost

%Final pertubations from the output X
D2=X.F2;
D1=X.F1;
D0=X.F0;

disp('FINAL Dist')
norm([D2-A2 D1-A1 D0-A0], 'fro')

function f = cost(X)
%Computes the cost function given an element X of the product manifold
    
   F2=X.F2;
   F1=X.F1;
   F0=X.F0;
   
   R = (F2 * V) * L^2 + (F1 * V) * L + (F0 * V);
      
   f=norm(R,'fro')^2;
end

function g = grad(X)
%Computes the euclidian gradient given an element X of the product manifold

   F2=X.F2;
   F1=X.F1;
   F0=X.F0;
      
   W2=V*L^2;
   W1=V*L; 
   W0=V;
   
   Matrix=F2*W2 + F1*W1 + F0*W0;
   
   g.F0=sparse_lr(F0, 2*Matrix, W0);
   g.F1=sparse_lr(F1, 2*Matrix, W1); 
   g.F2=sparse_lr(F2, 2*Matrix, W2);
   
   %In this way I usually compare the eucliedian gradient with the FD one
%    gradfun = approxgradientFD(problem);
%    temp1 = gradfun(X);
%    temp2 = problem.M.proj(X,g);
%    nn=norm(temp1.F2 - temp2.F2,'f')/ norm(temp2.F2,'fro')
%     
%     keyboard
end

end



