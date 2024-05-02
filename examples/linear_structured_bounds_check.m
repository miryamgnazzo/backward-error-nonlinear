%
% This script checks the structured bound that we have on randomly
% generated problems and compare this bound with the unstructured one
% (obtained with the Kathri Rao product).
%

ntests = 1000;
epsilon = 1e-3;
V = 0;
n = 64;

while cond(V) > 1e10
    F = cell(1, 5);

    F{1} = eye(n);
    F{2} = randn(n) / 10; F{2} = F{2} * F{2}';
    F{3} = randn(n); F{3} = -F{3} * F{3}';
    F{4} = randn(n); F{4} = F{4} * F{4}';
    F{5} = randn(n); F{5} = F{5} * F{5}';

    %Impose a certain sparsity pattern
    Sp2 = randi([0,1],n,n);
    Sp3 = randi([0,1],n,n);
    Sp4 = randi([0,1],n,n);
    Sp5 = randi([0,1],n,n);
    
    F{2} = F{2}.*Sp2;
    F{3} = F{3}.*Sp3;
    F{4} = F{4}.*Sp4;
    F{5} = F{5}.*Sp5;

    [V, L] = be_newton(F, @f, -1 : 1);

end

 %Construction of the matrix P for the linear structure
    P1=sparsity_struct(F{1});
    P2=sparsity_struct(F{2});
    P3=sparsity_struct(F{3});
    P4=sparsity_struct(F{4});
    P5=sparsity_struct(F{5});

    P=blkdiag(P1,P2,P3,P4,P5);

nrm = zeros(1, ntests);
be = zeros(1, ntests);
bnd = zeros(2, ntests);

for s = 1 : ntests
    Vt = V; Vt = Vt + randn(size(Vt)) * diag(epsilon * randn(1, size(Vt, 2)));
    Lt = L; Lt = Lt + diag(epsilon * randn(1, size(Lt, 1)));

    Ft = be_perturb(F, epsilon * exp(randn));
    
    Ft{1} =  (trace(Ft{1}) / n).*eye(n);
    Ft{2} = Ft{2}.*Sp2;
    Ft{3} = Ft{3}.*Sp3;
    Ft{4} = Ft{4}.*Sp4;
    Ft{5} = Ft{5}.*Sp5;

    R = be_residual(Ft, @f, V, L);
    nrm(s) = norm(R, 'fro');
    D = be_linear_structured(Ft, @f, V, L, P);

    be(s) = be_norm(D);

    %Unstructured bound 
    bnd(1, s) = be_unstructured_bound(3, Ft, @f, V, L);

    %Structured bound (linear structures)
    bnd(2, s) = be_linear_structured_bound(Ft, @f, V, L, P);
    
end

% Sort them
[nrm, I] = sort(nrm);
be = be(I);
bnd = bnd(:, I);

figure;
loglog(nrm, be, 'r*'); hold on;
plot(nrm, bnd(1, :), 'k--');
plot(nrm, bnd(2, :), 'b--');

 function [fv, fvp] = f(x)
     fv = [ x^2, x, 1, expm(-x), expm(-2*x)];
    if nargout > 1
        fvp = [2*x, 1, 0, -expm(-x), -2*expm(-2*x)];
    end
 end

 function P=sparsity_struct(A)
    %it reconstructs the matrix P associated with the matrix A
    %consider as structure the sparsity pattern induced by A
    
    [m,n]=size(A);
    [row,column,~]=find(A);
    
    s=length(row);
    
    P=zeros(m*n,s);
    
    for i = 1 : s
       
        index=(column(i)-1)*m + row(i);
    
        P(index,i) = 1; 
    end
end
