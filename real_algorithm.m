function [D0,D1,D2,e,t,infotable] = real_algorithm(A0, A1, A2, maxiter, timemax, Vhat, L,k)
%k is the rank of the matrix A1

n = length(A1);

mu=0; %10^-1;


%Manifold for the optimization
%problem.M = euclideancomplexfactory(n, 2*n-1);
%problem.M = spherecomplexfactory(n, m);
elements = struct();
elements.F2 = euclideanfactory(n,n); %coefficient of degree 2
%elements.F1 = euclideanfactory(n,n); %coefficient of degree 1
elements.F1=fixedrankembeddedfactory(n, n, k); %coefficient of degree 1
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

% %Starting point is [A0,A1,A2]
X0=M.rand();
X0.F2=A2;
[U1,S1,V1]=svd(A1);
X0.F1.U=U1(:,1:k);
X0.F1.S=(S1(1:k,1:k));
X0.F1.V=V1(:,1:k);
X0.F0=A0;

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

Y=M.randvec(X0);

%hessfd=getHessianFD(problem,X0,Y);
%keyboard

%Commenting the instruction for the gradient/Hessian we use the FD
%approximation in Manopt
%Uncommenting them we use the exact form provided in egrad or ehess

options.maxiter = maxiter;
options.maxtime = timemax;
options.tolgradnorm = 1e-12;
%options.Delta_bar = 4.47214*1e-0;
%options.Delta0 = options.Delta_bar/8;
options.debug=0;
options.rho_regularization = 1e3;

warning('on', 'manopt:getHessian:approx');

%START of the optimization (usually trustregions)
[X, xcost, info, options] = trustregions(problem,X0, options);
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
DD1=X.F1;
D0=X.F0;

%From structure to matrix (for the low rank coefficient)
D1=(DD1.U)*(DD1.S)*(DD1.V)';

disp('FINAL Dist')
norm([D2-A2 D1-A1 D0-A0],'fro')


function f = cost(X)
%Computes the cost function given an element X of the product manifold
    
   F2=X.F2;
   F1=X.F1;
   
   U=F1.U;
   V=F1.V;
   S=F1.S;
  
   F0=X.F0;
   
   FF=U*S*V'; %From structure to matrix for the low rank term
   
   F=[F2 FF F0];
   W=[Vhat*L^2; Vhat*L; Vhat]; 
   
      
 f=norm(F*W,'fro')^2 + mu*(norm(F-A,'fro')^2);

end

function g = egrad(X)
%Computes the euclidian gradient given an element X of the product manifold

   F2=X.F2;
   F1=X.F1;
  
   U=F1.U;
   V=F1.V;
   S=F1.S;
   
   FF1=U*S*V'; %From structure to matrix for the low rank term

   F0=X.F0;
      
   W2=Vhat*L^2; 
   W1=Vhat*L; 
   W0=Vhat;
   
   F=[F2 FF1 F0];
   W=[W2; W1; W0];
   
   Matrix=F*W;
   
   g.F0=2*Matrix*W0';
   g.F1=2*Matrix*W1'; 
   g.F2=2*Matrix*W2';
   
   g.F2=g.F2 + 2*mu*(F2-A2);
   g.F1=g.F1 + 2*mu*(FF1-A1);
   g.F0=g.F0 + 2*mu*(F0-A0);
   
   %In this way I usually compare the eucliedian gradient with the FD one
%    gradfun = approxgradientFD(problem);
%    temp1 = gradfun(X);
%    temp2 = problem.M.proj(X,g);
%        
%     nn=norm(temp1.F2 - temp2.F2,'f')/ norm(temp2.F2,'fro')
%     
%     keyboard
 end

end