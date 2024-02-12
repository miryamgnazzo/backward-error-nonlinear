%Script di prova per Backward error - Riemannian Optimization
%Consider a matrix polynomial \lambda^2 A_2 + \lambda A_1+ A_0

rng(2)
n=10; %size of the matrix polynomial
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

%Consider a perurbation of the exact eigenpairs -- only 2 eigenpairs here
%[VV,LL] = polyeig(Pert0,Pert1,Pert2); 

%[nn,I]=sort(abs(imag(LL)));
%numero di autovalori-autovettori approx
s=4;

V=VV(:,1:index(s));
L=(LL(1:index(s),1:index(s)));

[FF0,FF1,FF2]=unstruct_error(A0,A1,A2,V,L);


norm([FF0,FF1,FF2],'fro')
norm((A2+FF2)*V*L^2+ (A1+FF1)*V*L + (A0+FF0)*V,'fro')
keyboard

%taking real version (because fixedrankembeddedfactory is for real matrices)
%V=real(V);
%L=real(L);

%For Manopt
maxiter = 5000; %maximum number of iterations
timemax = 100; %maximum running time for manopt

%D0,D1,D2 are the minimizers of the functional
%[D0,D1,D2,e,t,infotable] = real_algorithm(A0, A1, A2, maxiter, timemax, V, L, k);


[D0,D1,D2]=struct_err(A0,A1,A2,V,L,k);


%check
norm((D2)*V*L^2+ (D1)*V*L + (D0)*V,'fro')

%norm([Pert2-A2 Pert1-A1 Pert0-A0],'fro')
norm([D2-A2 D1-A1 D0-A0],'fro')
