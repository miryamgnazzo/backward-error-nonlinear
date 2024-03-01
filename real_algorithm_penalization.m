function [D0,D1,D2,e,t,infotable] = real_algorithm_penalization(A0, A1, A2, I0, I1, I2, maxiter, timemax, Vhat, L, k, mu)
%k is the rank of the matrix A1

n = length(A1); 

%Manifold for the optimization
%problem.M = euclideancomplexfactory(n, 2*n-1);
%problem.M = spherecomplexfactory(n, m);
elements = struct();
elements.F2 = euclideanfactory(n,n); %coefficient of degree 2
%elements.F1 = euclideanfactory(n,n); %coefficient of degree 1
elements.F1 = fixedrankembeddedfactory(n, n, k); %coefficient of degree 1
elements.F0 = euclideanfactory(n,n); %coefficient of degree 0
M = productmanifold(elements); 

%Comment: we need a product manifold for optimize over the coefficients of
%the matrix polynomial. (F2,F1,F0) \in M: F2 and F0 are real matrices
%without additonal structures, F1 is a matrix of rank= k and real

problem.M=M;
problem.cost = @cost;

%Strating point is random
% X0=M.rand();
% X0.F2=randn(n);
% [U1,S1,V1]=svd(randn(n));
% X0.F1.U=U1(:,1:k);
% X0.F1.S=(S1(1:k,1:k));
% X0.F1.V=V1(:,1:k);
% X0.F0=randn(n);

% % %Starting point is [A0,A1,A2]
% X0=M.rand();
% X0.F2=A2;
% [U1,S1,V1]=svd(A1);
% X0.F1.U=U1(:,1:k);
% X0.F1.S=(S1(1:k,1:k));
% X0.F1.V=V1(:,1:k);
% X0.F0=A0;


%Strating point is given by [I0, I1, I2]
X0 = M.rand();
X0.F2 = I2;
[U1,S1,V1] = svd(I1);
X0.F1.U = U1(:,1:k);
X0.F1.S = (S1(1:k,1:k));
X0.F1.V = V1(:,1:k);
X0.F0 = I0;


%starting point from unstruct
% [FF0,FF1,FF2]=unstruct_error(A0,A1,A2,Vhat,L);
% X0=M.rand();
% X0.F2=A2+FF2;
% [U1,S1,V1]=svd(A1+FF1);
% X0.F1.U=U1(:,1:k);
% X0.F1.S=(S1(1:k,1:k));
% X0.F1.V=V1(:,1:k);
% X0.F0=A0+FF0;

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

function f = cost(X)
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

end