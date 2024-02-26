%Script unstructured

rng(2)
n=10; %size of the matrix polynomial
k=2; %rank of the coefficient A1

A2=randn(n);
%A1=randn(n);
U1=randn(n,k);
V1=rand(k,n);
A1=U1*V1;
A0=randn(n);

epsilon=10^-3;
Delta2=epsilon*randn(n);
Delta0=epsilon*randn(n);
Delta1_U1=epsilon*randn(n,k);
Delta1_V1=epsilon*randn(k,n);

Pert2=A2+Delta2;
Pert0=A0+Delta0;
Pert1=(U1+Delta1_U1)*(V1+Delta1_V1);

[VV,LL,index]=approx(Pert0,Pert1,Pert2);

s=2;

V=VV(:,1:index(s));
L=(LL(1:index(s),1:index(s)));

[D0,D1,D2]=unstruct_error(A0,A1,A2,V,L);
norm((A2+D2)*V*L^2+ (A1+D1)*V*L + (A0+D0)*V,'fro')
