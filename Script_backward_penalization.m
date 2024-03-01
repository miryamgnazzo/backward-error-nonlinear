% Script for the numerical approximation of the Backward error
% Considering a penalization term \mu 

%Script di prova per Backward error - Riemannian Optimization
%Consider a matrix polynomial \lambda^2 A_2 + \lambda A_1+ A_0

rng(4)
n=100; %size of the matrix polynomial
k=2; %rank of the coefficient A1

A2=randn(n);
%A1=randn(n);
U1=randn(n,k);
V1=rand(k,n);
A1=U1*V1;
A0=randn(n);

epsilon=10^-1;
Delta2=epsilon*randn(n);
Delta0=epsilon*randn(n);
Delta1_U1=epsilon*randn(n,k);
Delta1_V1=epsilon*randn(k,n);

Pert2=A2+Delta2;
Pert0=A0+Delta0;
Pert1=(U1+Delta1_U1)*(V1+Delta1_V1);

[VV,LL,index]=approx(Pert0,Pert1,Pert2);

%sottospazio invariante al posto di metter i complessi per gli autovettori
%mackey de teran dopico quasicanonical form for quadratic matirx
%polynomials

%numero di autovalori-autovettori approx
s=4;

V=VV(:,1:index(s));
L=(LL(1:index(s),1:index(s)));

[FF0,FF1,FF2]=unstruct_error(A0,A1,A2,V,L);

norm([FF0,FF1,FF2],'fro')
norm((A2+FF2)*V*L^2+ (A1+FF1)*V*L + (A0+FF0)*V,'fro')

%For Manopt
maxiter = 100; %maximum number of iterations
timemax = 10; %maximum running time for manopt

mu=10.^(1:10);
%mu=10.^(-(1:9));

I0 = A0;
I1 = A1;
I2 = A2;
 
% I0 = randn(n);
% I1 = randn(n,k)*randn(k,n);
% I2 = randn(n);

f1 = mu(1)*norm(A2*V*L^2+ A1*V*L + A0*V,'fro')^2 + norm([A0,A1,A2],'fro')^2;
%f1 = norm(A2*V*L^2+ A1*V*L + A0*V,'fro')^2 + mu(1)*norm([A0,A1,A2],'fro')^2;

for i = 1: length(mu)
    %I_i starting points
    [D0,D1,D2,e,t,infotable] = real_algorithm_penalization(A0, A1, A2, I0, I1, I2, maxiter, timemax, V, L, k, mu(i));

 f2 = e(end);

%     if (abs(f1-f2)<1e-3*mu(i))
%         
%         fprintf('Number iteration %d with mu %e \n', i, mu(i))
%         break;
% 
%     else

%  keyboard

        I0 = D0;
        I1 = D1;
        I2 = D2;
        
        f1 = f2;

%    end

    
end

%check
fprintf('Check correct %e \n',norm((A2+D2)*V*L^2+ (A1+D1)*V*L + (A0+D0)*V,'fro'))
%norm((A2+D2)*V*L^2+ (A1+D1)*V*L + (A0+D0)*V,'fro')


fprintf('Structured %e > Unstructured %e \n',norm([D2 D1 D0],'fro'), norm([FF2 FF1 FF0],'fro'))
norm([D2 D1 D0],'fro')